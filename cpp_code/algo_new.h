#ifndef FASTIFDS_ALGO_NEW_H
#define FASTIFDS_ALGO_NEW_H

#include "template.h"
#include "IfdsInstance.h"
#include "LCA.h"



ull FIRST_32BIT = (1ULL << 32) - 1;
int HashModG;

struct HashG {
	size_t operator() (const pii& p) const {
		return p.fs * 1LL * HashModG + p.sc;
	}
};


class algo_new {
public:


	const int START_VERTEX = 0;
	const int EXIT_VERTEX = 1;
	const int CALL_VERTEX = 2;
	const int RETURN_SITE_VERTEX = 3;

	void stopClock(clock_t Time, string s = "") {

		cout << "[" << s << "]" << " time taken = " <<
		     double(clock() - Time) / CLOCKS_PER_SEC << " s" << endl;
	}


	IfdsInstance instance;
	// Copied from the instance
	int n_G, n_H, m_G, m_H, n_GExp, m_GExp;
	vector<pair<int, int>> edgeListG, edgeListH;
	vector<vector<int>> G, H, GExp;
	vector<vector<vector<int>>> flow;
	vector<int> vertexTypeG, procOf;
	vector<int> D;
	vector<vector<vector<int>>> TWD;
	vector<int> TWD_root;
	vector<vector<int>> TWD_par;
	vector<int> par_tdH;
	vector<vector<int>> nodeGExp;
	vector<pair<int, int>> revMapGExp;
	vvi nodesPerProc;


	// new stuff
	V<V<V<V<int>>>> callerEdges, retSiteEdges;  // as in the other algos
	int n_HExp;
	LCA tdLCAHExp;
	vector<LCA> TWDLCA;
	vector<vector<int>> GHat, GExpRev;
	unordered_set<pair<int, int>, HashG>  edgeListGHat;
	vector<unordered_set<ull>>  edgeListGHatPerProc;
	vector<int> s;  // for i in [0, n_H), s[i] in [0, n_G) is the start node of proc_i in G (NOT GExp)
	vector<int> e;  // for i in [0, n_H), e[i] in [0, n_G) is the exit node of proc_i in G (NOT GExp)
	vector<vector<vector<int>>> TWD_adj;    // TWD[i] for i in [0, n_H) is graph corresponding to T_i
	vector<vector<int>> order;
	vector<unordered_set<pair<int, int>, HashG>> RLocalPerProc; // RLocalPerProc[p] holds edges (U, V) in GHat_p
	vector<vector<int>> twd_Depth;   // twd_Depth[p][j] = depth of bag j in [0, n_{B_p}) in treedecomp. of p.
	vector<vector<int>> delta;   // delta[p][b] = sum of sizes of bag b and all of its ancestors
	vector<unordered_set<ull>> RAncPerProc; // RAncPerProc[p] holds edges (U, V) in GHat_p
	vector<int> rb;     // rb[u] for u in [0, n_G) is a bag label in [0, |TWD[procOf[u]]|), which is
	// the root bag of u. i.e. the lowest-depth bag that has u.
	vector<vector<int>> bagsContainingNode; // bagsContainingNode[u] is all bags in [0, |TWD[procOf[u]]|) that contain u.[
	vector<vector<vector<long long>>> RAncBit;    // RAncBit[u][d] is a bit vector of delta[rb[u]]*(D[procOf[u]]+1) bits
	// for ancestor-bags of rb[u] in order of INCREASING depth: we add |D+1| bits for each node in the bag
	// if a bit is 1, means there is an edge from (u, d) to the node in the corresponding bit position
	vector<vector<vector<long long>>> RAncBitRev;
	// Same as RAncBit, but if a bit is 1, means there is an edge from the node in the corresponding bit position TO (u, d)

	vvi callsGExp; // callsGExp[U] = array of start-nodes V that are reachable from U by first taking a same-proc path in
					// in GExp then taking a direct edge to V. Both U and V are nodes of GExp.

	vi maxBagSz;

	// H comes in
	vector<vector<int>> HExp, HExpRev;
	vector<vector<int>> nodeHExp;
	vector<pair<int, int>> revMapHExp;
	unordered_set<pair<int, int>, HashG> edgeListHExp;
	vector<int> par_tdHExp; // we build a treedepth dec. of HExp out of our treedepth dec. of H
	vector<vector<int>> tdHExpAdj;
	int tdHExpRoot;
	vector<gp_hash_table<int, null_type>> Sin, Sout, anc_hash;  // anc_hash is just differing way of storing anc
	vector<int> tdd_Depth;
	vector<vector<bool>> f, f_R;    // f[U][j] for j in [0, tdd_Depth[U]) is 1 if U's ancestor of depth j is in f according
	// to f's definition in the paper.
	vector<vector<int>> anc;    // anc[U] for U in HExp is the ancestors of U in treedepth dec. of HExp (not including U)

	clock_t timeDiff;

	algo_new(IfdsInstance _instance) : instance(_instance) {

		copyData();


		n_HExp = 0;

		for (int i = 0; i < n_H; ++i)
			n_HExp += D[i];

		HashModG = n_GExp;

		doAuxiliaryStuff();

		preprocess();
	}





	void preprocess() {

		/*
		 * Step 0: finding callsGExp
		 */

		clock_t Time = clock();

		buildCallsGExp();


		stopClock(Time, "Step 0: finding callsGExp");

		/*
		 * Step 3: finding LCA for treewidth decompositions
		 */

		Time = clock();

		TWDLCA.assign(n_H, LCA());

		for (int i = 0; i < n_H; ++i)
			TWDLCA[i] = LCA(TWD_par[i]);


		stopClock(Time, "Step 3: finding LCA for treewidth decompositions");

		/*
		 * Step 4: Computing GHat. Algorithm 1 in the ESOP paper
		 * Takes a few seconds
		 */
		Time = clock();

		buildGHat();

		stopClock(Time, "Step 4: Computing GHat");

		/*
		 * Step 5: Local preprocessing
		 */
		Time = clock();


		computeTWD_adj();


		timeDiff = 0;

		order.assign(n_H, vi());
		RLocalPerProc.assign(n_H, unordered_set<pair<int, int>, HashG>());
		edgeListGHatPerProc.assign(n_H, unordered_set<ull>());

		for (auto& pr : edgeListGHat) {
			assert(procOf[revMapGExp[pr.fs].fs] == procOf[revMapGExp[pr.sc].fs]);
			edgeListGHatPerProc[procOf[revMapGExp[pr.fs].fs]].insert(((ull (pr.fs))  << 32) | pr.sc);
		}

		// also compute twd_Depth, delta, rb, bagsContainingNode on the fly
		twd_Depth.assign(n_H, vi());
		delta.assign(n_H, vi());
		rb.assign(n_G, -1);
		bagsContainingNode.assign(n_G, vi());

		// easy part, needn't be parallelized
		for (int p = 0; p < n_H; ++p) {
			twd_Depth[p].assign(sz(TWD[p]), 0);
			delta[p].assign(sz(TWD[p]), 0);
			dfsRLocal(TWD_adj[p], TWD_root[p], -1, p, 0, 0);
		}


		validateTWD();  //  auxiliary, just to check that TWD is valid


		// Calls to solve doesn't use GHat, so we'll clear it and use results of solve to re-build its
		// up-to-date version
		for (int U = 0; U < n_GExp; ++U)
			GHat[U].clear();

		maxBagSz.assign(n_H, 0);


		for (int p = 0; p < n_H; ++p) {
			for (auto& bag : TWD[p]) {
				maxBagSz[p] = max(maxBagSz[p], sz(bag));
			}
		}


		int threadNum = 35;


		double tot = 0;

		for (int p = 0; p < n_H; ++p) {
			tot += double(sz(nodesPerProc[p])) * maxBagSz[p] * maxBagSz[p] * maxBagSz[p] * D[p] * D[p] * D[p];
		}

		db(tot);



		vvi parts;
		V<double> partsSz;

		parts.pb(vi());
		partsSz.pb(0);


		double cur = 0;

		for (int p = 0; p < n_H; ++p) {
			parts.back().pb(p);
			cur += double(sz(nodesPerProc[p])) * maxBagSz[p] * maxBagSz[p] * maxBagSz[p] * D[p] * D[p] * D[p];
			partsSz.back() = cur;
//			db(cur, tot);
			if (cur * threadNum >= tot)
				parts.pb(vi()), partsSz.pb(0), cur = 0;
		}

		if (parts.back().empty())
			parts.pop_back(), partsSz.pop_back();

		db(sz(parts), partsSz);

		omp_set_num_threads(sz(parts));
		#pragma omp parallel
		{

			int nthreads = omp_get_num_threads();
			int processor_id = omp_get_thread_num();
			cout << 0 << flush;

			for (int p : parts[processor_id]) {
				solve(p, 0);
				for (ull UV : edgeListGHatPerProc[p])
					GHat[UV >> 32].pb(UV & FIRST_32BIT);
			}
			cout << 1 << flush;

		}
		cout << endl;



		// IMP: from now on, only use edgeListGHatPerProc not edgeListGHat

		int RLocalSz = 0;
		for (int p = 0; p < n_H; ++p)
			RLocalSz += sz(RLocalPerProc[p]);

		db(RLocalSz);

		stopClock(Time, "Step 5: Local preprocessing in same-bag nodes");

		/*
		 * Step 6: Ancestor preprocessing
		 */
		Time = clock();

		int mxBagDepth = 0;
		for (int p = 0; p < n_H; ++p) {
			for (int b = 0; b < sz(TWD[p]); ++b) {
				mxBagDepth = max(mxBagDepth, twd_Depth[p][b]);
			}
		}
		db(mxBagDepth);

		RAncPerProc.assign(n_H, unordered_set<ull>());




		omp_set_num_threads(sz(parts));
		#pragma omp parallel
		{

			int nthreads = omp_get_num_threads();
			int processor_id = omp_get_thread_num();
			cout << 0 << flush;

			for (int p : parts[processor_id]) {
				computeAncReachability(p);
			}
			cout << 1 << flush;
		}
		cout << endl;

		int RAncSz = 0;

		for (int p = 0; p < n_H; ++p) {
			RAncSz += sz(RAncPerProc[p]);
		}

		db(RAncSz);

		stopClock(Time, "Step 5.5: Ancestor preprocessing before RAncBit");

		Time = clock();

		buildRAncBit();


		stopClock(Time, "Step 6: Ancestor preprocessing");

		/*
		 * Step 7: Computing HExp
		 */
		Time = clock();

		buildHExp();
//		buildHExpDumb();


		stopClock(Time, "Step 7: Computing HExp");
		/*
		 * Step 8: Building par_tdHExp and its LCA
		 */
		Time = clock();


		buildPar_tdHExp();



		tdHExpAdj.assign(n_HExp, vi());

		for (int i = 0; i < n_HExp; ++i)
			if (par_tdHExp[i] != -1)
				tdHExpAdj[i].pb(par_tdHExp[i]), tdHExpAdj[par_tdHExp[i]].pb(i);


		// we'll need HExpRev

		HExpRev.assign(n_HExp, vi());
		for (pair<int, int> edge : edgeListHExp)
			HExpRev[edge.sc].pb(edge.fs);


		anc_hash.assign(n_HExp, gp_hash_table<int, null_type>());
		anc.assign(n_HExp, vi());

		buildAnc(tdHExpRoot, vector<int>());
		int mxAncSz = 0;
		for (int U = 0; U < n_HExp; ++U) {
			assert(U == tdHExpRoot || !anc[U].empty());
			mxAncSz = max(mxAncSz, sz(anc[U]) + 1);
		}
		db(mxAncSz);

		validateTDD();



		tdLCAHExp = LCA(par_tdHExp);

		stopClock(Time, "Step 8: Building par_tdHExp and its LCA");
		/*
		 * Step 9: Building S
		 */
		Time = clock();
		buildS();
		stopClock(Time, "Step 9: Building S");
		/*
		 * Step 10: Building F
		 */
		Time = clock();
		buildF();
		stopClock(Time, "Step 10: Building F");
	}

	// all the book-keeping data-structures that don't rely on any step outside this function
	void doAuxiliaryStuff() {


		GExpRev.assign(n_GExp, vi());
		for (int U = 0; U < n_GExp; ++U) {
			for (int V : GExp[U])
				GExpRev[V].pb(U);
		}


		callerEdges.assign(n_H, VVV<int>());
		retSiteEdges.assign(n_H, VVV<int>());

		for (int c = 0; c < n_G; ++c) {
			if (vertexTypeG[c] == CALL_VERTEX) {
				int r = -1;
				assert(sz(G[c]) == 2);
				for (auto& v : G[c]) {
					if (vertexTypeG[v] == RETURN_SITE_VERTEX)
						assert(r == -1), r = v;
				}
				assert(r != -1);
				assert(procOf[c] == procOf[r]);

				int p = -1;  // called procedure
				for (auto& v : G[c]) {
					if (vertexTypeG[v] == START_VERTEX)
						assert(p == -1), p = procOf[v];
				}
				assert(p != -1);


				callerEdges[p].pb(VV<int>(D[p] + 1));

				for (auto& C : nodeGExp[c]) {
					for (auto& S : GExp[C]) {
						int s_p = revMapGExp[S].fs, d = revMapGExp[S].sc;
						if (vertexTypeG[s_p] == START_VERTEX) {
							callerEdges[p].back()[d].pb(C);
						}
					}
				}

				retSiteEdges[p].pb(VV<int>(D[p] + 1));

				for (auto& R : nodeGExp[r]) {
					for (auto& E : GExpRev[R]) {
						int e_p = revMapGExp[E].fs, d = revMapGExp[E].sc;
						if (vertexTypeG[e_p] == EXIT_VERTEX) {
							retSiteEdges[p].back()[d].pb(R);
						}
					}
				}

			}
		}



		for (int p = 0; p < n_H; ++p)
			assert(s[p] != -1 && e[p] != -1);



	}

	void buildCallsGExp() {
		callsGExp.assign(n_GExp, vi());

		vi vis(n_GExp, 0);
		int visCnt = 0;

		for (int C = 0; C < n_GExp; ++C) {
			int c = revMapGExp[C].fs;
			if (vertexTypeG[c] == CALL_VERTEX) {
				for (auto& S : GExp[C]) {
				    int s = revMapGExp[S].fs;
					if (vertexTypeG[s] == START_VERTEX) {
						++visCnt;

						function<void(int)> dfs = [&](int U) {
							callsGExp[U].pb(S);
							vis[U] = visCnt;

							for (auto& V : GExpRev[U])
								if (procOf[revMapGExp[V].fs] == procOf[c] && vis[V] != visCnt)
									dfs(V);

						};


						dfs(C);
					}
				}
			}
		}

		// making the elements of callsGExp[U] distinct
		for (int U = 0; U < n_GExp; ++U) {
			si s;
			for (auto& elem : callsGExp[U]) {
			    s.insert(elem);
			}
			callsGExp[U].clear();
			for (auto& elem : s) {
				callsGExp[U].pb(elem);
			}
		}

		int totSz = 0;
		int mxSz = 0;

		for (int U = 0; U < n_GExp; ++U) {
			totSz += sz(callsGExp[U]);
			mxSz = max(mxSz, sz(callsGExp[U]));
		}

	}

	void buildGHat() {

		GHat.assign(n_GExp, vi());


		// Lines 1-3:
		unordered_set<ull> Q, S, EDash;
		Q.reserve(1 << 25);
		S.reserve(1 << 25);
		EDash.reserve(1 << 25);


		vector<vector<int>> EDashAdj(n_GExp);

		{
			int cnt = 0;
			for (int U = 0; U < n_GExp; ++U) {
				for (int V : GExp[U])
					++cnt;
			}
			db(cnt);
		}
		for (int U = 0; U < n_GExp; ++U) {
			for (int V : GExp[U])
				Q.insert((((ull) U) << 32) | V);
		}

		db(sz(Q));

		int cnt = 0;
		for (int i = 0; i < 100000; ++i) {
			cnt += present(Q, rand());
		}

		int iterationCnt = 0;

		while (!Q.empty()) {

//			db(sz(Q), sz(EDash), sz(Q) + sz(EDash));
			// Lines 5-10
			ull UV = *Q.begin();
			++iterationCnt;
			if (iterationCnt % 1000000 == 0)
				db(iterationCnt);
			Q.erase(UV);
			int U = UV >> 32, V = UV & FIRST_32BIT;
			int u = revMapGExp[U].fs, d1 = revMapGExp[U].sc;
			int v = revMapGExp[V].fs, d2 = revMapGExp[V].sc;

			bool isInterProcedural = (vertexTypeG[u] == CALL_VERTEX && vertexTypeG[v] == START_VERTEX)
			                         || (vertexTypeG[u] == EXIT_VERTEX && vertexTypeG[v] == RETURN_SITE_VERTEX);
			if (isInterProcedural)
				continue;

			int p = procOf[u];

			EDash.insert(UV);
			EDashAdj[U].pb(V);

			// Lines 11-13

			for (int d3 = 0; d3 <= D[p]; ++d3) {
				int S = nodeGExp[s[p]][d3];
				if (present(EDash, (((ull) S) << 32) | U) || S == U) {
					if (!present(EDash, (((ull) S) << 32) | V) && !present(Q, (((ull) S) << 32) | V))
						Q.insert((((ull) S) << 32) | V);
				}
			}

			// Lines 14-17

			assert((u == s[p]) == (vertexTypeG[u] == START_VERTEX));
			if (u == s[p]) {
				for (int W : EDashAdj[V]) {
					if (!present(EDash, (((ull) U) << 32) | W) && !present(Q, (((ull) U) << 32) | W))
						Q.insert((((ull) U) << 32) | W);
				}
			}

			// Lines 18-23
			// Lie 19-20 are involving EBar not EDash (i.e. GExp)
			if (u == s[p] && v == e[p]) {
				for (int i = 0; i < sz(callerEdges[p]); ++i) {
					for (auto& C : callerEdges[p][i][d1]) {
						for (auto& R : retSiteEdges[p][i][d2]) {
							assert(procOf[revMapGExp[C].fs] == procOf[revMapGExp[R].fs]);
							if (!present(EDash, (((ull) C) << 32) | R) && !present(Q, (((ull) C) << 32) | R))
								S.insert((((ull) C) << 32) | R), Q.insert((((ull) C) << 32) | R);
						}
					}
				}
			}

		}

		// Lines 24-28


		for (int U = 0; U < n_GExp; ++U) {
			for (int V : GExp[U]) {
				int u = revMapGExp[U].fs;
				int v = revMapGExp[V].fs;
				// TODO: check if disallowing self-loops in GHat will make any difference
				if (procOf[u] == procOf[v])
					edgeListGHat.insert(mp(U, V));
			}
//			edgeListGHat.insert(mp(U, U));
		}

		for (ull UV : S)
			edgeListGHat.insert(mp(UV >> 32, UV & FIRST_32BIT));

		for (pair<int, int> edge : edgeListGHat) {
			int U = edge.fs, V = edge.sc;
			GHat[U].pb(V);
		}

		db(sz(edgeListGHat));

	}


	void computeAncReachability(int p) {
		int elemCntInF = 0;
		int n_B = sz(TWD[p]);
		V<V<V<V<unsigned long long>>>> F(n_B), FDash(n_B);
		// F[b][i][d][dep]: u := i'th node in bag b, corresponds to F(u, d, b, dep)
		// 0 <= dep <= depth of bag b.
		// Here we do the bit-packing as described in the paper
		for (int b = 0; b < n_B; ++b) {
			int bagSz = sz(TWD[p][b]);
			F[b].assign(bagSz, VV<unsigned long long>());
			FDash[b].assign(bagSz, VV<unsigned long long>());
			for (int i = 0; i < bagSz; ++i) {
				int Sz = delta[p][b] * (D[p] + 1);
				F[b][i].assign(D[p] + 1, V<unsigned long long>((Sz + 63) / 64));
				FDash[b][i].assign(D[p] + 1, V<unsigned long long>((Sz + 63) / 64));
				elemCntInF += (D[p] + 1) * delta[p][b] * (D[p] + 1);
			}
		}
		vi order;
		for (int i = 0; i < n_B; ++i) {
			order.pb(i);
		}
		sort(order.begin(), order.end(), [&](int b1, int b2) {
			return twd_Depth[p][b1] < twd_Depth[p][b2];
		});

		for (auto& b : order) {
			int b_p = TWD_par[p][b];
			if (b_p != -1) {
				unordered_set<int> b_pSet;
				unordered_map<int, int> b_pNodeToIdx;
				for (int j = 0; j < sz(TWD[p][b_p]); ++j) {
					int u = TWD[p][b_p][j];
					b_pSet.insert(u);
					b_pNodeToIdx[u] = j;
				}

				for (int i = 0; i < sz(TWD[p][b]); ++i) {
					int u = TWD[p][b][i];
					if (present(b_pSet, u)) {
						for (int d = 0; d <= D[p]; ++d) {
							int j = b_pNodeToIdx[u];
							int Sz_p = delta[p][b_p] * (D[p] + 1);
							int shift = sz(TWD[p][b]) * (D[p] + 1);
							for (int k = 0; k < Sz_p; ++k) {
								// set (k+shift)'th bit in F[b][i][d]
								// to be k'th bit in F[b_p][j][d], same for FDash
								int A_k = k / 64, B_k = k % 64;
								int A_ksh = (k + shift) / 64, B_ksh = (k + shift) % 64;

								F[b][i][d][A_ksh] |= ((F[b_p][j][d][A_k] >> B_k) & 1ULL) << B_ksh;
								FDash[b][i][d][A_ksh] |= ((FDash[b_p][j][d][A_k] >> B_k) & 1ULL) << B_ksh;
							}
						}
					}
				}
			}


			for (int i = 0; i < sz(TWD[p][b]); ++i) {
				int u = TWD[p][b][i];
				for (int d1 = 0; d1 <= D[p]; ++d1) {
					int U = nodeGExp[u][d1];
					for (int j = 0; j < sz(TWD[p][b]); ++j) {
						int v = TWD[p][b][j];
						for (int d2 = 0; d2 <= D[p]; ++d2) {
							int V = nodeGExp[v][d2];
							if (present(RLocalPerProc[p], mp(U, V))) {
								int bitPosF = j * (D[p] + 1) + d2, bitPosFDash = i * (D[p] + 1) + d1;
								F[b][i][d1][bitPosF / 64] |= 1ULL << (bitPosF % 64);
								FDash[b][j][d2][bitPosFDash / 64] |= 1ULL << (bitPosFDash % 64);

								int Sz = sz(F[b][i][d1]);
								assert(Sz == sz(FDash[b][j][d2]));
								// Do ORing
								int BagBlockSz = sz(TWD[p][b]) * (D[p] + 1);

								if (BagBlockSz % 64) {
									int div = BagBlockSz / 64, mod = BagBlockSz % 64;
									F[b][i][d1][div] |= (F[b][j][d2][div] >> (mod)) << mod;
									FDash[b][j][d2][div] |= (FDash[b][i][d1][div] >> (mod)) << mod;
								}
								for (int k = (BagBlockSz + 63) / 64; k < Sz; ++k) {
									F[b][i][d1][k] |= F[b][j][d2][k];
									FDash[b][j][d2][k] |= FDash[b][i][d1][k];
								}
							}
						}
					}
				}
			}

		}

//		db(elemCntInF);


		for (int b = 0; b < n_B; ++b) {
			int bagSz = sz(TWD[p][b]);
			for (int i = 0; i < bagSz; ++i) {
				int u = TWD[p][b][i];
				for (int d = 0; d <= D[p]; ++d) {
					int U = nodeGExp[u][d];
					int curBagIdx = b;
					int startIdx = 0;
					for (int j = 0; j <= twd_Depth[p][b]; ++j) {

						int Sz = sz(TWD[p][curBagIdx]) * (D[p] + 1);

						for (int k = 0; k < Sz; ++k) {
							// (k+startIdx)'th bit being true in F[b][i][d] means that (TWD[p][curBagIdx][k / D], k % DStar) is
							// stored in F[b][i][d][twd_Depth[p][curBagIdx]], same for FDash

							int kSt = k + startIdx;
							int A_kSt = kSt / 64, B_kSt = kSt % 64;

							int V = nodeGExp[TWD[p][curBagIdx][k / (D[p] + 1)]][k % (D[p] + 1)];

							if ((F[b][i][d][A_kSt] >> B_kSt) & 1ULL) {
								RAncPerProc[p].insert((((ull) U) << 32) | V);
							}


							if ((FDash[b][i][d][A_kSt] >> B_kSt) & 1ULL) {

								RAncPerProc[p].insert((((ull) V) << 32) | U);
							}

						}

						startIdx += Sz;
						curBagIdx = TWD_par[p][curBagIdx];
					}
					assert(curBagIdx == -1);
					assert(startIdx + 64 > sz(F[b][i][d]) * 64);
					assert(startIdx <= sz(F[b][i][d]) * 64);
				}
			}
		}


	}

	void solve(int p, int idx) {

		// Lines 5-9

		int bagLabel = order[p][idx];
		vi bag = TWD[p][bagLabel];
		int bagSz = sz(bag);
		vector<int> vertexSet;

		assert(sz(bag) <= 12);

		for (int u : bag)
			for (int d = 0; d <= D[procOf[u]]; ++d) {
				int U = nodeGExp[u][d];
				vertexSet.pb(U);
			}

		vector<pair<int, int>> edges = computeAllPairsReachability(vertexSet, p);

		for (pair<int, int> edge : edges) {
			edgeListGHatPerProc[p].insert(((ull (edge.fs))  << 32) | edge.sc);
			RLocalPerProc[p].insert(edge);
		}

		// Lines 10-14

		if (idx + 1 < sz(order[p])) {
			solve(p, idx + 1);

			vector<pair<int, int>> edges = computeAllPairsReachability(vertexSet, p);

			for (pair<int, int> edge : edges) {
				edgeListGHatPerProc[p].insert(((ull (edge.fs))  << 32) | edge.sc);
				RLocalPerProc[p].insert(edge);
			}
		}
	}

	void dfsRLocal(vvi& adj, int u, int parent, int p, int dep, int currDelta) {

		for (auto& v : TWD[p][u]) {
			if (rb[v] == -1)
				rb[v] = u;
			bagsContainingNode[v].pb(u);
		}

		twd_Depth[p][u] = dep;
		currDelta += sz(TWD[p][u]);
		delta[p][u] = currDelta;
		for (int v : adj[u])
			if (v != parent)
				dfsRLocal(adj, v, u, p, dep + 1, currDelta);
		assert(p < sz(order));
		order[p].pb(u);
	}

	void computeTWD_adj() {
		TWD_adj.assign(n_H, vvi());

		for (int i = 0; i < n_H; ++i) {

			assert(sz(TWD[i]) == sz(TWD_par[i]));
			int n_T_i = sz(TWD[i]); // number of bags in treew. dec. of G_i
			TWD_adj[i].assign(n_T_i, vi());

			for (int j = 0; j < n_T_i; ++j)
				if (TWD_par[i][j] != -1) {
					TWD_adj[i][j].pb(TWD_par[i][j]);
					TWD_adj[i][TWD_par[i][j]].pb(j);
				}
		}
	}

	vector<pair<int, int>> computeAllPairsReachability(vector<int> vertexSet, int p) {
		int N = sz(vertexSet);

//		db(N);
		vector<vector<int>> adj(N);
		vector<pair<int, int>> ret;
		vector<bool> vis(N);

		int e = 0;

		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				if (present(edgeListGHatPerProc[p], ((ull (vertexSet[i]))  << 32) | vertexSet[j]))
					adj[i].pb(j), ++e;

//		db(e, N * N);


		function<void(int, int)> dfsLocal = [&](int src, int u) {

			vis[u] = true;
			ret.pb(mp(vertexSet[src], vertexSet[u]));

			for (auto& v : adj[u]) {
				if (!vis[v])
					dfsLocal(src, v);
			}
		};

		for (int s = 0; s < N; ++s)
			vis.assign(N, 0), dfsLocal(s, s);

		return ret;

	}


	void buildRAncBit() {

		// Build RAncBit
		RAncBit.assign(n_G, VV<ll>());
		RAncBitRev.assign(n_G, VV<ll>());

		vvi bagIds(n_G); // bagIds[u] = [rb[u], par(rb[u]), par(par(rb[u])), .. , root]

		// nodeToIdxPerProc[p][i][u] = k means that in TWD of p, in i'th bag, node u appears and appears at TWD[p][i][k]
		VV<unordered_map<int, int>> nodeToIdxPerProc(n_H);

//		#pragma omp parallel for
		for (int p = 0; p < n_H; ++p) {
			int n_B = sz(TWD[p]);
			nodeToIdxPerProc[p].assign(n_B, unordered_map<int, int>());
			for (int i = 0; i < sz(TWD[p]); ++i) {
				for (int j = 0; j < sz(TWD[p][i]); ++j) {
					nodeToIdxPerProc[p][i][TWD[p][i][j]] = j;
				}
			}
		}


//		#pragma omp parallel for
		for (int u = 0; u < n_G; ++u) {
			int p = procOf[u];
			RAncBit[u].assign(D[p] + 1, V<ll>());
			RAncBitRev[u].assign(D[p] + 1, V<ll>());

			for (int d = 0; d <= D[p]; ++d) {
				int U = nodeGExp[u][d];
				int bitSz = (D[p] + 1) * delta[p][rb[u]];
				RAncBit[u][d].assign((bitSz + 63) / 64, 0);
				RAncBitRev[u][d].assign((bitSz + 63) / 64, 0);
			}


			vi bags;
			int curBag = rb[u];
			do {
				bags.pb(curBag);
				curBag = TWD_par[p][curBag];
			} while (curBag != -1);
			reverse(bags.begin(), bags.end());

			bagIds[u] = bags;
		}

		int tot = 0;


//		#pragma omp parallel for
		for (int p = 0; p < n_H; ++p) {
			for (auto& UV : RAncPerProc[p]) {
				int U = UV >> 32, V = UV & FIRST_32BIT;
				int u = revMapGExp[U].fs, d1 = revMapGExp[U].sc;
				int v = revMapGExp[V].fs, d2 = revMapGExp[V].sc;
				tot += twd_Depth[p][rb[u]] + twd_Depth[p][rb[u]];

				int curIdx = 0;

				for (auto& b_id : bagIds[u]) {
					if (present(nodeToIdxPerProc[p][b_id], v)) {
						int bitIdx = curIdx + nodeToIdxPerProc[p][b_id][v] * (D[p] + 1) + d2;
						RAncBit[u][d1][bitIdx / 64] |= (1ULL << (bitIdx % 64));
					}
					curIdx += (D[p] + 1) * sz(TWD[p][b_id]);
				}

				curIdx = 0;

				for (auto& b_id : bagIds[v]) {
					if (present(nodeToIdxPerProc[p][b_id], u)) {
						int bitIdx = curIdx + nodeToIdxPerProc[p][b_id][u] * (D[p] + 1) + d1;
						RAncBitRev[v][d2][bitIdx / 64] |= (1ULL << (bitIdx % 64));
					}
					curIdx += (D[p] + 1) * sz(TWD[p][b_id]);
				}
			}
		}
	}




	bool querySameContext(int U, int V) {
		int u = revMapGExp[U].fs, v = revMapGExp[V].fs;
		if (procOf[u] != procOf[v])
			return 0;
		int p = procOf[u];
		int bu = rb[u], bv = rb[v], b = TWDLCA[p].query(bu, bv);

		bool ret = false;

		for (auto& w : TWD[p][b]) {
			for (int d = 0; d <= D[p]; ++d) {
				int W = nodeGExp[w][d];
				// IMP: & vs &&, the & to make things fairer
				if (present(RAncPerProc[p], (((ull) U) << 32) | W) & present(RAncPerProc[p], (((ull) W) << 32) | V))
					ret = true;
			}
		}
		return ret;
	}

	bool querySameContextFast(int U, int V) {
		int u = revMapGExp[U].fs, v = revMapGExp[V].fs;
		if (procOf[u] != procOf[v])
			return 0;
		int d1 = revMapGExp[U].sc, d2 = revMapGExp[V].sc;
		int p = procOf[u];
		int bu = rb[u], bv = rb[v], b = TWDLCA[p].query(bu, bv), dep = twd_Depth[p][b];

		// Look at bag dep'th bag in RAncBit, RAncBitRev
		int idxSt = (delta[p][b] - sz(TWD[p][b])) * (D[p] + 1), idxEn = delta[p][b] * (D[p] + 1) - 1;
		// for idx in [idxSt, idxEn]: if idx'th bit in RAncBit[u][d1] = RAncBitRev[v][d2] = 1: return true

		int divSt = idxSt / 64, modSt = idxSt % 64;
		int divEn = idxEn / 64, modEn = idxEn % 64;

		if (divSt == divEn)
			return (((RAncBit[u][d1][divSt] & RAncBitRev[v][d2][divSt]) >> modSt) & ((1ULL << (modEn - modSt + 1)) - 1)) > 0;

		if (modSt) {
			if ((RAncBit[u][d1][divSt] & RAncBitRev[v][d2][divSt]) >> modSt)
				return true;
			modSt = 0, ++divSt;
		}
		if (modEn != 63) {
			if ((RAncBit[u][d1][divEn] & RAncBitRev[v][d2][divEn]) & ((1ULL << (modEn + 1)) - 1))
				return true;
			modEn = 63, --divEn;
		}

		bool ret = false;

		for (int i = divSt; i <= divEn; ++i) {
			if (RAncBit[u][d1][i] & RAncBitRev[v][d2][i])
				ret = true;
		}

		return ret;
	}

	bool queryUisStart(int U, int V) {
		// handles only the case where u is a start nodes
		int u = revMapGExp[U].fs, v = revMapGExp[V].fs;
		int d1 = revMapGExp[U].sc, d2 = revMapGExp[V].sc;
		assert(vertexTypeG[u] == START_VERTEX);

		bool ret = 0;

		if (vertexTypeG[v] == START_VERTEX) {

			if (procOf[u] == procOf[v])
				return querySameContextFast(U, V);

			int UU = nodeHExp[procOf[u]][d1], VV = nodeHExp[procOf[v]][d2];

			int WW = tdLCAHExp.query(UU, VV);
			if (WW == UU)
				return f_R[VV][tdd_Depth[UU]];
			if (WW == VV)
				return f[UU][tdd_Depth[VV]];

			for (int dep = 0; dep <= tdd_Depth[WW]; ++dep) {
				if (f[UU][dep] && f_R[VV][dep]) {
					ret = true;
				}
			}
		} else {
			for (int d3 = 0; d3 <= D[procOf[v]]; ++d3) {
				int W = nodeGExp[s[procOf[v]]][d3];
				if (querySameContextFast(W, V) & queryUisStart(U, W)) {
					ret = true;
				}
			}
		}
		return ret;

	}

	bool query(int U, int V) {
		int u = revMapGExp[U].fs, v = revMapGExp[V].fs;
		int d1 = revMapGExp[U].sc, d2 = revMapGExp[V].sc;

		if (vertexTypeG[u] == START_VERTEX)
			return queryUisStart(U, V);

		bool ret = false;
		for (auto& S : callsGExp[U]) {
		    if (queryUisStart(S, V))
				ret = true;
		}

		return ret;
	}

	void buildHExpDumb() { // build a maximal HExp, just to measure performance
		HExp.assign(n_HExp, vi());
		for (int p = 0; p < n_H; ++p) {
			for (auto& pp : H[p]) {
				for (int d1 = 0; d1 <= D[p]; ++d1) {
					for (int d2 = 0; d2 <= D[pp]; ++d2) {
						int U = nodeHExp[p][d1], V = nodeHExp[pp][d2];
						HExp[U].pb(V);
						edgeListHExp.insert(mp(U, V));
					}
				}
			}
		}
	}
	void buildHExp() {

		// Building nodeHExp, revMapHexp as specified in the comment in buildHExp

		n_HExp = 0;
		nodeHExp.assign(n_H, vi());

		for (int p = 0; p < n_H; ++p)
			for (int d = 0; d <= D[p]; ++d)
				nodeHExp[p].pb(n_HExp++), revMapHExp.pb(mp(p, d));

		/*
		 * s[p] is a node in G which is the start node of procedure p.
		 *
		 * HExp: Exploded graph of H. Defined as in the paper. It has n_HExp := sum_{p in H} (D[p]+1) nodes
		 * Similar to GExp in IfdsInstance.h, we'll have vvi nodeHExp(n_H) where nodeHExp[p] is a vector of size
		 * D[p]+1 with the label of nodes in the 'exploded' p (w.r.t. H NOT G).
		 * Their order is consistent with the order flow[e] for all e in [0, m_G) where s[p] appears.
		 *
		 * We'll label the nodes of HExp from 0 to n_HExp := (sum_{p in H} (D[p]+1))-1, and we'll also store a reverse
		 * mapping vector<pii> revMapHExp(n_GExp) where:
		 *             revMapHExp[x] = (p, d) (for x in [0, n_HExp), p in [0, n_H), d in [0, D[p]])
		 *
		 *
		 * We have an edge from (p_1, d_1) to (p_2, d_2) iff it's possible to start from (s[p_1], d_1) in GHat,
		 * take some same-context-path, then reach a call node directly going out to (s[p_2], d_2).
		 *
		 * for each node U := (u, d) in GExp where u is of type CALL_VERTEX:
		 *      edgeDestVector = [] // Destination of edges that will be added to edgeListHExp, will have O(1) size
		 *      for each of its neighbours in GExp (v, d`):
		 *          if (v is of type START_VERTEX):
		 *              edgeDestVector.pb((procOf[v], d`))
		 *      p = procOf[u]
		 *      for d_s in [0, D[p]]
		 *          V := (p, d_s)
		 *          ask a querySameContext on GHat for (s[p], d_s)->(u, d)
		 *          if it's true:
		 *              for each W in edgeDestVector:
		 *                  add (V, W) to HExp
		 *
		 */

		for (int U = 0; U < n_GExp; ++U) {
			int u = revMapGExp[U].fs, d1 = revMapGExp[U].sc;
			if (vertexTypeG[u] != CALL_VERTEX)
				continue;
			vector<int> startFromCall;
			for (int V : GExp[U]) {
				int v = revMapGExp[V].fs, d2 = revMapGExp[V].sc;
				if (vertexTypeG[v] == START_VERTEX)
					startFromCall.pb(nodeHExp[procOf[v]][d2]);
			}
			int p = procOf[u];
			for (int d_s = 0; d_s <= D[p]; ++d_s) {
				int V = nodeHExp[p][d_s];
				if (querySameContextFast(nodeGExp[s[p]][d_s], nodeGExp[u][d1])) {

					// IMP: check if removing the condition will make a difference, it shouldn't
					//  this is because same-proc edges in HExp are accounted for by same-context queries
					for (int W : startFromCall)
						edgeListHExp.insert(mp(V, W));
				}
			}
		}

		HExp.assign(n_HExp, vi());

		for (pair<int, int> edge : edgeListHExp)
			HExp[edge.fs].pb(edge.sc);

	}

	void buildPar_tdHExp() {

		par_tdHExp.assign(n_HExp, -1);
		tdHExpRoot = -1;
		for (int i = 0; i < n_H; ++i)
			if (par_tdH[i] == -1) {
				assert(tdHExpRoot == -1);
				tdHExpRoot = nodeHExp[i][0];
			}

		assert(tdHExpRoot != -1);

		for (int i = 0; i < n_H; ++i) {
			for (int j = 1; j <= D[i]; ++j)
				par_tdHExp[nodeHExp[i][j]] = nodeHExp[i][j - 1];
			if (nodeHExp[i][0] != tdHExpRoot)
				par_tdHExp[nodeHExp[i][0]] = nodeHExp[par_tdH[i]][D[par_tdH[i]]];
		}
	}

	void buildAnc(int U, vector<int> ancestors, int parent = -1) {
		for (int ancestor : ancestors)
			anc_hash[U].insert(ancestor);
		anc[U] = ancestors;
		ancestors.pb(U);

		for (int V : tdHExpAdj[U])
			if (V != parent)
				buildAnc(V, ancestors, U);
	}


	// initially U = subTreeRoot, do dfsRLocal from U in adj, traversing a node V only iff present(anc_hash[V], subTreeRoot)
	// add all the nodes you visit to S (which means we can also use it as a vis array)
	void dfsS(int U, int subTreeRoot, vector<vector<int>>& adj, vector<gp_hash_table<int, null_type>>& S) {
		S[subTreeRoot].insert(U);
		for (int V : adj[U])
			if (!present(S[subTreeRoot], V) && present(anc_hash[V], subTreeRoot))
				dfsS(V, subTreeRoot, adj, S);
	}


	void buildS() {


		tdd_Depth.assign(n_HExp, 0);
		for (int U = 0; U < n_HExp; ++U)
			tdd_Depth[U] = sz(anc[U]);


		Sin.assign(n_HExp, gp_hash_table<int, null_type>());
		Sout.assign(n_HExp, gp_hash_table<int, null_type>());

		// computing S_in[U], S_out[U]
		for (int U = 0; U < n_HExp; ++U) {
			Sin[U].insert(U);
			Sout[U].insert(U);
			// Find nodes that U can reach
			dfsS(U, U, HExp, Sout);
			// Find nodes that can reach U
			dfsS(U, U, HExpRev, Sin);
		}
	}

	void buildF() {
		/*
		 * Computing f and f_R, this part can be parallelized
		 */


		f.assign(n_HExp, vector<bool>());
		f_R.assign(n_HExp, vector<bool>());

		int maxTreedepth_Depth = 0;

		for (int i = 0; i < n_HExp; ++i) {
			maxTreedepth_Depth = max(maxTreedepth_Depth, tdd_Depth[i]);
		}

		ll tot = 0;
		int threadNum = 20;

		for (int U = 0; U < n_HExp; ++U) {
			tot += tdd_Depth[U] * 1LL * tdd_Depth[U];
		}

		vvi parts;

		parts.pb(vi());

		ll cur = 0;

		for (int U = 0; U < n_HExp; ++U) {
			cur += tdd_Depth[U] * 1LL * tdd_Depth[U];
			if (cur * threadNum >= tot)
				parts.pb(vi()), cur = tdd_Depth[U] * 1LL * tdd_Depth[U];
			parts.back().pb(U);
		}

		db(sz(parts));

		omp_set_num_threads(sz(parts));
#pragma omp parallel
		{

			int nthreads = omp_get_num_threads();
			int processor_id = omp_get_thread_num();
			cout << 0 << flush;

			for (int U : parts[processor_id]) {
				f[U].assign(tdd_Depth[U] + 1, 0);
				f_R[U].assign(tdd_Depth[U] + 1, 0);

				for (int i = 0; i < tdd_Depth[U]; ++i) {
					int V = anc[U][i];
					for (int j = 0; j <= i; ++j) {
						int W = anc[U][j];
						f[U][i] = f[U][i] | (present(Sin[W], U) && present(Sout[W], V));
						f_R[U][i] = f_R[U][i] | (present(Sout[W], U) && present(Sin[W], V));
					}
				}
			}
			cout << 1 << flush;
		}

		cout << endl;
	}

	void validateTWD() {
		for (int u = 0; u < n_G; ++u) {
			sort(bagsContainingNode[u].begin(), bagsContainingNode[u].end(), [&](int b1, int b2) {
				return twd_Depth[procOf[u]][b1] < twd_Depth[procOf[u]][b2];
			});
		}



		// check that bagsContainingNode[u] is connected and non-empty

		for (int u = 0; u < n_G; ++u) {
			si prevBags;
			assert(!bagsContainingNode[u].empty());
			assert(rb[u] != -1);
			assert(bagsContainingNode[u][0] == rb[u]);
			prevBags.insert(rb[u]);
			int p = procOf[u];

			for (auto& b : bagsContainingNode[u]) if (b != rb[u]) {
					assert(present(prevBags, TWD_par[p][b]));
					prevBags.insert(b);
				}
		}


		// check that each edge appears in a bag

		for (int u = 0; u < n_G; ++u) {
			int p = procOf[u];
			si canHave; // set of all neighbours u can have given that TWD is valid
			for (auto& b : bagsContainingNode[u]) {
				for (auto& v : TWD[p][b]) {
					canHave.insert(v);
				}
			}
			for (auto& v : G[u]) {
				if (procOf[v] == procOf[u])
					assert(present(canHave, v));
			}
		}
		int mxBagSz = 0;

		for (int p = 0; p < n_H; ++p) {
			for (auto& elem : TWD[p]) {
			    mxBagSz = max(mxBagSz, sz(elem));
			}
		}

		cout << "TWD is valid! with max bag size = " << mxBagSz << endl;
	}

	void validateTDD() { // checking that the constructed treedepth decomp. of HExp is valid
		int cnt = 0;
		function<void(int, int)> dfsHExpT = [&](int u, int par) {
			++cnt;
			for (auto& v : tdHExpAdj[u]) {
				if (v != par)
					dfsHExpT(v, u);
			}
		};

		dfsHExpT(tdHExpRoot, -1);
		assert(cnt == n_HExp);


		for (int U = 0; U < n_HExp; ++U) {
			for (auto& V : HExp[U])
				assert(present(anc_hash[U], V) || present(anc_hash[V], U) || U == V);
		}

		cout << "TDD is valid!" << endl;
	}

	void copyData() {
		n_G = instance.n_G;
		n_H = instance.n_H;
		n_GExp = instance.n_GExp;
		m_G = instance.m_G;
		m_H = instance.m_H;
		m_GExp = instance.m_GExp;
		D = instance.D;
		edgeListG = instance.edgeListG;
		edgeListH = instance.edgeListH;
		G = instance.G;
		H = instance.H;
		GExp = instance.GExp;
		flow = instance.flow;
		vertexTypeG = instance.vertexTypeG;
		procOf = instance.procOf;
		TWD = instance.TWD;
		TWD_root = instance.TWD_root;
		TWD_par = instance.TWD_par;
		par_tdH = instance.par_tdH;
		nodeGExp = instance.nodeGExp;
		revMapGExp = instance.revMapGExp;
		s = instance.s;
		e = instance.e;
		nodesPerProc = instance.nodesPerProc;
	}



};

#endif //FASTIFDS_ALGO_NEW_H
