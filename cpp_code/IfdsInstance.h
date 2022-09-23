#ifndef FASTIFDS_IFDSINSTANCE_H
#define FASTIFDS_IFDSINSTANCE_H

#include "template.h"

class IfdsInstance {
public:
	const int START_VERTEX = 0;
	const int EXIT_VERTEX = 1;
	const int CALL_VERTEX = 2;
	const int RETURN_SITE_VERTEX = 3;
	const int BAD_PROC_THRESHOLD = 50000;

	/*
	 * Each instance corresponds to a codebase, and in our case each benchmark from the DaCapo benchmark suite maps to
	 * one instance.
	 *
	 * Note: nothing is built artificially here
	 *
	 * An instance consists of:
	 *      - Directed graph G of n_G nodes labeled with [0..n_G) and m_G edges labeled with [0..m_G)
	 *          - Describes the CFGs of all n_H procedures, which are labeled with [0, n_H). G_i for i in [0, n_H)
	 *          corresponds to the CFG of i'th procedure. For G_0, its vertex set is labeled with [0, n_{G_0}),
	 *          and the vertex set of G_1 is labeled with [n_{G_0}, n_{G_0}+n_{G_1}) and so on..
	 *              - Each procedure has:
	 *                  - Unique start vertex marked with vertexTypeG[u] = START_VERTEX
	 *                  - Unique exit vertex marked with vertexTypeG[u] = EXIT_VERTEX
	 *              - For inter-procedural calls, each call is represented by:
	 *                  - Unique call vertex marked with vertexTypeG[u] = CALL_VERTEX
	 *                  - Unique return-site vertex marked with vertexTypeG[u] = RETURN_SITE_VERTEX
	 *                  - The call vertex has exactly two outgoing edges: call-to-return-site and call-to-start
	 *                  - The return-site vertex has exactly two incoming edges: call-to-return-site and exit-to-return-site
	 *
	 *      - Directed graph H of n_H nodes labeled with [0..n_H) and m_H edges labeled with [0..m_G)
	 *          - Describes the call graph of G: there is an edge (i, j) if somewhere in G_i we have a call to G_j
	 *
	 *      - For G and H, we are given their edge list, and an edge's label is its index in the list.
	 *
	 *      - Treewidth decomposition of n_H CFGs G_i:
	 *              - TWD[i] for i in [0, n_H) is an array of vectors of length n_{TW_i}
	 *              - Bags of G_i are labeled with [0, n_{TW_i}), TWD[i][j] j'th bag in the treewidth decomposition of G_i
	 *              - TWD_root[i] for i in [0, n_H) is in [0, n_{TW_i}): the root of treewidth decomposition of G_i
	 *              - TWD_par[i][j] for i in [0, n_H) and j in [0, n_{TW_j}) is parent of j'th bag of twd of G_i
	 *              and is -1 only in one root node.
	 *
	 *
	 *      - Treedepth decomposition of H: an array par_tdH[i] for i in [0, n_H) where par_tdH[i] = -1 if i is a root
	 *      (there is one unique root, the dec. is connected), otherwise par_tdH[i] is the parent of i in the treedepth dec.
	 *
	 *      - An array flow[i] for i in [0, m_G) is a vector of pairs describing the flow function associated with i'th
	 *          edge in G. A pair (a, b) (for a, b in [0, DStar)) in M[i] corresponds to a UP-DOWN edge in the succinct
	 *          representations of the flow function. Its format is according to FastIFDS.java
	 *
	 *      - D := D_0, D_1, .., D_{n_H-1}, denoting the domain size of each function
	 *
	 * From these, we can easily construct the following and give it to the solvers as part of its input:
	 *
	 *      - GExp: Exploded graph of G. Defined as in the paper. We'll use the following format to label nodes of GExp:
	 *      we'll have vvi nodeGExp(n_G) where nodeGExp[u] is a vector of size D[procOf[u]]+1 with the label of nodes
	 *      in the 'exploded' u. Their order is consistent with the order flow[e] for all e in [0, m_G) where u appears.
	 *
	 *      We'll label the nodes of GExp from 0 to n_GExp := (sum_{u in G} (D[procOf[u]]+1))-1, and we'll also store a reverse
	 *      mapping vector<pii> revMapGExp(n_GExp) where:
	 *              revMapGExp[x] = (u, d) (for x in [0, n_GExp), u in [0, n_G), d in [0, D[procOf[u]]])
	 */
	int n_G, n_H, m_G, m_H, n_GExp, m_GExp;
	vector<pair<int, int>> edgeListG, edgeListH;
	vector<vector<int>> G, H, GExp;
	vector<vector<vector<int>>> flow;
	vector<int> vertexTypeG, procOf; // for i in [0, n_G), procOf[i] is in [0, n_H)
	vector<int> D;
	vector<vector<vector<int>>> TWD;
	vector<int> TWD_root;
	vector<vector<int>> TWD_par;
	vector<int> par_tdH;
	vector<vector<int>> nodeGExp;
	vector<pair<int, int>> revMapGExp;
	vector<int> s;  // for i in [0, n_H), s[i] in [0, n_G) is the start node of proc_i in G (NOT GExp)
	vector<int> e;  // for i in [0, n_H), e[i] in [0, n_G) is the exit node of proc_i in G (NOT GExp)
	vvi nodesPerProc;
	vector<bool> badProc;   // badProc[p] = true if CFGSz[p] * D[p] > BAD_PROC_THRESHOLD
	int badProcCnt;
	double badProcPercentage;
	vvi callsGExp;  // same as in algo_new
	vector<vector<int>> GExpRev;

	IfdsInstance() {}



	IfdsInstance(int nG, int nH, int mG, int mH, vector<int>& d, const vector<pair<int, int>> &edgeListG,
	             const vector<pair<int, int>> &edgeListH, const vector<vector<vector<int>>> &flow,
	             const vector<int> &vertexTypeG, const vector<int> &procOf, const vector<vector<vector<int>>> &twd,
	             const vector<int> &twdRoot, const vector<vector<int>> &twdPar, const vector<int> &parTdH) : n_G(nG),
				 n_H(nH), m_G(mG), m_H(mH), D(d), edgeListG(edgeListG), edgeListH(edgeListH), flow(flow),
				 vertexTypeG(vertexTypeG), procOf(procOf), TWD(twd), TWD_root(twdRoot), TWD_par(twdPar), par_tdH(parTdH) {
		buildNodesPerProc();
		buildG();
		buildH();
		buildStAndEx();
		buildGExp();

		buildBadProc();
		eliminateBadProc();


		buildNodesPerProc();
		buildG();
		buildH();
		buildStAndEx();
		buildGExp();

		db(n_GExp, m_GExp);
	}



	void buildNodesPerProc() {
		nodesPerProc.assign(n_H, vi());
		for (int u = 0; u < n_G; ++u) {
			nodesPerProc[procOf[u]].pb(u);
		}
	}


	void eliminateBadProc() {

		/*
		 * int nG, int nH, int m_G, int mH, vector<int>& d, const vector<pair<int, int>> &edgeListG,
	             const vector<pair<int, int>> &edgeListH, const vector<vector<vector<int>>> &flow,
	             const vector<int> &vertexTypeG, const vector<int> &procOf, const vector<vector<vector<int>>> &twd,
	             const vector<int> &twdRoot, const vector<vector<int>> &twdPar, const vector<int> &parTdH
		 */
		int Tn_G, Tn_H;     // Done
		int Tm_G, Tm_H;     // Done
		vi TD;              // Done
		vii TedgeListG;     // Done
		vii TedgeListH;     // Done
		V<vvi> Tflow;       // Done
		vi TvertexTypeG;    // Done
		vi TprocOf;         // Done
		VVV<int> TTWD;      // Done
		vi TTWD_root;       // Done
		vvi TTWD_par;       // Done
		vi Tpar_tdH;        // Done


		vi ICFGMapping(n_G, -1);
		vi methodMapping(n_H, -1);

		Tn_G = Tn_H = 0;


		for (int i = 0; i < n_H; ++i) {
			if (!badProc[i])
				methodMapping[i] = Tn_H++;
		}

		for (int i = 0; i < n_G; ++i)
			if (!badProc[procOf[i]])
				ICFGMapping[i] = Tn_G++;


		TprocOf.assign(Tn_G, -1);

		for (int i = 0; i < n_G; ++i) if (ICFGMapping[i] != -1)
			TprocOf[ICFGMapping[i]] = methodMapping[procOf[i]];

		for (int i = 0; i < Tn_G; ++i) {
			assert(TprocOf[i] != -1);
			assert(TprocOf[i] < Tn_H);
		}


		int eId = 0;

		for (auto& pr : edgeListG) {
			int u = pr.fs, v = pr.sc;
			int sat = (ICFGMapping[u] != -1) + (ICFGMapping[v] != -1);
			if (sat == 2) {
				TedgeListG.pb(mp(ICFGMapping[u], ICFGMapping[v]));

				vvi newFlow = flow[eId];
				assert(sz(newFlow) == D[procOf[u]]  + 1);



				// if (u, v) is a call-to-return-site and the called function is bad, set flow to identity
				if (procOf[u] == procOf[v] && vertexTypeG[u] == CALL_VERTEX && vertexTypeG[v] != START_VERTEX) {
					assert(vertexTypeG[v] == RETURN_SITE_VERTEX);

					int start = -1;
					for (auto& w : G[u]) {
					    if (w != v)
							assert(start == -1), start = w;
					}
					assert(start != -1);
					if (badProc[procOf[start]]) {
						for (int i = 0; i <= D[procOf[u]]; ++i) {
							newFlow[i].clear();
							newFlow[i].pb(i);
						}
						vertexTypeG[u] = vertexTypeG[v] = -1;
					}
				}
				Tflow.pb(newFlow);
			}
			++eId;
		}


		for (auto& pr : edgeListH) {
			int p1 = pr.fs, p2 = pr.sc;
			int sat = (methodMapping[p1] != -1) + (methodMapping[p2] != -1);
			if (sat == 2)
				TedgeListH.pb(mp(methodMapping[p1], methodMapping[p2]));
		}



		TvertexTypeG.assign(Tn_G, -2);

		for (int i = 0; i < n_G; ++i) if (ICFGMapping[i] != -1)
			TvertexTypeG[ICFGMapping[i]] = vertexTypeG[i];


		Tm_G = sz(TedgeListG), Tm_H = sz(TedgeListH);

		for (int i = 0; i < Tn_G; ++i) {
			assert(TvertexTypeG[i] != -2);
		}


		TD.assign(Tn_H, -1);

		for (int i = 0; i < n_H; ++i) if (methodMapping[i] != -1)
			TD[methodMapping[i]] = D[i];

		for (int i = 0; i < Tn_H; ++i)
			assert(TD[i] != -1);


		TTWD.assign(Tn_H, vvi());
		TTWD_par.assign(Tn_H, vi());
		TTWD_root.assign(Tn_H, 0);

		for (int i = 0; i < n_H; ++i)  if (methodMapping[i] != -1) {

			// Has id's from G, need to update them
			TTWD[methodMapping[i]] = TWD[i];
			for (int j = 0; j < sz(TTWD[methodMapping[i]]); ++j) {
				for (int k = 0; k < sz(TTWD[methodMapping[i]][j]); ++k) {
					assert(ICFGMapping[TTWD[methodMapping[i]][j][k]] != -1);
					TTWD[methodMapping[i]][j][k] = ICFGMapping[TTWD[methodMapping[i]][j][k]];
				}
			}

			// unchanged, since this concerns labeling of bags within a CFG, which is 0-based for all methods.
			TTWD_par[methodMapping[i]] = TWD_par[i];
			TTWD_root[methodMapping[i]] = TWD_root[i];

		}




		Tpar_tdH.assign(Tn_H, -2);
		int newRoot = -1;

		for (int p = 0; p < n_H; ++p) if (methodMapping[p] != -1) {
			assert(!badProc[p]);
			int cur = par_tdH[p];
			while (cur != -1 && badProc[cur])
				cur = par_tdH[cur];
			if (cur == -1)
				Tpar_tdH[methodMapping[p]] = -1, newRoot = methodMapping[p];
			else
				assert(methodMapping[cur] != -1), Tpar_tdH[methodMapping[p]] = methodMapping[cur];
		}

		assert(newRoot != -1);

		for (int i = 0; i < Tn_H; ++i) {
			assert(Tpar_tdH[i] != -2);
			if (Tpar_tdH[i] == -1 && i != newRoot)
				Tpar_tdH[i] = newRoot;
		}

		int minusOnes = 0;

		for (int i = 0; i < Tn_H; ++i) {
			if (Tpar_tdH[i] == -1)
				++minusOnes;
		}

		assert(minusOnes == 1);

		n_G = Tn_G, n_H = Tn_H;
		m_G = Tm_G, m_H = Tm_H;
		D = TD;
		edgeListG = TedgeListG;
		edgeListH = TedgeListH;
		flow = Tflow;
		vertexTypeG = TvertexTypeG;
		procOf = TprocOf;
		TWD = TTWD;
		TWD_par = TTWD_par;
		TWD_root = TTWD_root;
		par_tdH = Tpar_tdH;

	}

	void buildBadProc() {
		badProc.assign(n_H, 0);
		badProcCnt = 0;

		vi callsCnt(n_H);

		for (int u = 0; u < n_G; ++u) {
			if (vertexTypeG[u] == CALL_VERTEX)
				++callsCnt[procOf[u]];
		}

		for (int p = 0; p < n_H; ++p)
			if (sz(nodesPerProc[p]) * 1LL * D[p] > BAD_PROC_THRESHOLD) {
				badProc[p] = true, ++badProcCnt;
				db(D[p], sz(nodesPerProc[p]));
			}




		badProcPercentage = badProcCnt * 100.0 / n_H;

		for (int p = 0; p < n_H; ++p) {
			if (badProc[p]) {
				map<int, bool> vis;
				int cycles = 0;
				mii order;
				function<void(int)> dfs = [&](int u) {
					order[u] = sz(order);
					vis[u] = true;
					for (auto& v : G[u]) if (procOf[u] == procOf[v]) {
						if (!vis[v])
					        dfs(v);
						else
							++cycles;
					}
				};
//				dfs(s[p]);
//				db(sz(order));
//				db(cycles);
			}
		}

	}


	void buildStAndEx() {
		s.assign(n_H, -1);
		e.assign(n_H, -1);
		for (int u = 0; u < n_G; ++u) {
			if (vertexTypeG[u] == START_VERTEX) {
				assert(s[procOf[u]] == -1);
				s[procOf[u]] = u;
			}
			if (vertexTypeG[u] == EXIT_VERTEX) {
				assert(e[procOf[u]] == -1);
				e[procOf[u]] = u;
			}
		}
		for (int p = 0; p < n_H; ++p) {
			assert(s[p] != -1 && e[p] != -1);
		}
	}

	void buildG() {
		G.assign(n_G, vector<int>());

		for (auto &pr: edgeListG)
			G[pr.fs].pb(pr.sc);
	}

	void buildH() {
		H.assign(n_H, vector<int>());

		for (auto &pr: edgeListH)
			H[pr.fs].pb(pr.sc);
	}


	void buildGExp() {

		n_GExp = 0;
		m_GExp = 0;
		nodeGExp.assign(n_G, vi());
		revMapGExp.clear();

		for (int i = 0; i < n_G; ++i)
			for (int d = 0; d <= D[procOf[i]]; ++d)
				nodeGExp[i].pb(n_GExp++), revMapGExp.pb(mp(i, d));

		GExp.assign(n_GExp, vector<int>());

		for (int i = 0; i < m_G; ++i) {
			int u = edgeListG[i].fs, v = edgeListG[i].sc;
			assert(u != v);
			for (int d1 = 0; d1 <= D[procOf[u]]; ++d1)
				for (auto& d2 : flow[i][d1])
				    GExp[nodeGExp[u][d1]].pb(nodeGExp[v][d2]), ++m_GExp;
		}

	}

	void buildCallsGExp() {
		GExpRev.assign(n_GExp, vi());
		for (int U = 0; U < n_GExp; ++U) {
			for (int V : GExp[U])
				GExpRev[V].pb(U);
		}

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


	map<int, vi> generateQuerySet() {

		buildCallsGExp();

		map<int, vi> ret;

		// generate n_G random (U, V)'s but making sure callsGExp[U] is non-empty
		int Sz = 0;
		while (Sz < n_G) {
			int U = rand() % n_GExp;
			if (!callsGExp[U].empty())
				ret[U].pb(rand() % n_GExp), ++Sz;;
		}

		db(sz(ret) * 1.0 / n_G);

		return ret;

/*		int cnt = n_G / n_H;
		for (int p = 0; p < n_H; ++p) {
			// IMP
			vi Us({nodeGExp[s[p]][0], nodeGExp[s[p]][1]});

			for (auto& U : Us) {

				assert(vertexTypeG[revMapGExp[U].fs] == START_VERTEX);

				assert(vertexTypeG[s[p]] == START_VERTEX);


				for (int j = 0; j < cnt; ++j) {
					int V = rand() % n_GExp;
					ret[p][U].pb(V);
				}


				vi closeProcedures, nodesInCloseProcedures;
				{
					vi dist(n_H, 10000000);
					queue<int> q;
					q.push(p);
					dist[p] = 0;
					while (!q.empty()) {
						int u = q.front();
						closeProcedures.pb(u);
						if (sz(closeProcedures) == 20)
							break;
						q.pop();
						for (auto& v : H[u])
							if (dist[v] == 10000000)
								q.push(v), dist[v] = dist[u] + 1;
					}

					for (auto& pp : closeProcedures) {
						for (auto& u : nodesPerProc[pp]) {
							nodesInCloseProcedures.pb(u);
						}
					}
				}

				for (int j = 0; j < cnt; ++j) {
					int v = nodesInCloseProcedures[rand() % sz(nodesInCloseProcedures)];
					int d = rand() % (D[procOf[v]] + 1);
					int V = nodeGExp[v][d];
					ret[p][U].pb(V);
				}

				for (int j = 0; j < cnt; ++j) {
					int v = nodesPerProc[p][rand() % sz(nodesPerProc[p])];
					int d = rand() % (D[p] + 1);
					int V = nodeGExp[v][d];
					ret[p][U].pb(V);
				}

			}
		}

		return ret;*/
	}
};


#endif //FASTIFDS_IFDSINSTANCE_H
