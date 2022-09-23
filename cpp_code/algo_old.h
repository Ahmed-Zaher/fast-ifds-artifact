#ifndef FASTIFDS_ALGO_OLD_H
#define FASTIFDS_ALGO_OLD_H

#include "template.h"
#include "SimpleIfdsInstance.h"

int HashModG_alg_old;

struct HashG_alg_old {
	size_t operator() (const pii& p) const {
		return p.fs * HashModG_alg_old + p.sc;
	}
};


class algo_old {
public:
	void stopClock(clock_t Time, string s = "") {

		cout << "[" << s << "]" << " time taken = " <<
		     double(clock() - Time) / CLOCKS_PER_SEC << " s" << endl;
	}

	const int START_VERTEX = 0;
	const int EXIT_VERTEX = 1;
	const int CALL_VERTEX = 2;
	const int RETURN_SITE_VERTEX = 3;

	int n_G, n_H, n_GExp;
	vector<vector<int>> G, GExp;
	vector<int> vertexTypeG, procOf;
	vector<int> D;
	vector<vector<int>> nodeGExp;
	vector<pair<int, int>> revMapGExp;

	int mainNodeInGExp; // The U in a (U, V) query.
	vi s;               // s[p] for p in [0, n_H) is a node in G which is the start of proc p.
	vi returnSite;      // returnSite[u] for u in [0, n_G) and u is a CALL_VERTEX, is the corresponding
	// RETURN_SITE_VERTEX in G.

	// The rest is direct translation of the paper, line numbering is according to the original version
	// (we make the variables global instead of passing them all the time as in the conference-version)

	unordered_set<ull> PathEdge, WorkList, SummaryEdge;
	vvi SummaryEdgeAdj, PathEdgeAdjRev;
	vvi GExpRev;
	V<V<V<V<int>>>> callerEdges, retSiteEdges;

	ull FIRST_32BIT;

	V<si> X;  // X[u] for u in [0, n_G) is a set of elements taking values in [0, D[procOf[u]]], which
	// are the set of data facts that hold at node u

	clock_t timeForForwardTabulate;
	clock_t runTime;
	clock_t Time;




	algo_old(SimpleIfdsInstance* instance, int _mainNodeInGExp)
			: mainNodeInGExp(_mainNodeInGExp) {
		Time = clock();

		FIRST_32BIT = (1ULL << 32) - 1;

		{
			n_G = instance->n_G;
			n_H = instance->n_H;
			n_GExp = instance->n_GExp;
			D = instance->D;
			G = instance->G;
			GExp = instance->GExp;
			vertexTypeG = instance->vertexTypeG;
			procOf = instance->procOf;
			nodeGExp = instance->nodeGExp;
			revMapGExp = instance->revMapGExp;
			s = instance->s;
		}
		int u = revMapGExp[mainNodeInGExp].fs;
		int p = procOf[u];

		if (vertexTypeG[u] != START_VERTEX) {
			GExp[nodeGExp[s[p]][0]].clear();
			GExp[nodeGExp[s[p]][0]].pb(mainNodeInGExp);
			mainNodeInGExp = nodeGExp[s[p]][0];
		}

		HashModG_alg_old = n_GExp;
		buildReturnSite();

		SummaryEdgeAdj.assign(n_GExp, vi());
		PathEdgeAdjRev.assign(n_GExp, vi());
		GExpRev.assign(n_GExp, vi());

		for (int U = 0; U < n_GExp; ++U) {
			for (auto& V : GExp[U]) {
				GExpRev[V].pb(U);
			}
		}

		callerEdges.assign(n_H, VVV<int>());
		retSiteEdges.assign(n_H, VVV<int>());

		// Suppose for procedure p, it has calls labeled starting from 0. Each call has
		// a call-node and a return-site-node from the caller's procedure.
		// callerEdges[p][j][d] records the edges (c, d2) -> (s[p], d) for different d2 where
		// c is a call-node involved in the j'th call to procedure p.
		// retSiteEdges[p][j][d] records the edges (e[p], d) -> (r, d2) for different d2 where
		// r is a return-site-node involved in the j'th call to procedure p.
		for (int c = 0; c < n_G; ++c) {
			if (vertexTypeG[c] == CALL_VERTEX) {
				int r = returnSite[c];
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

//		db("Did the work before SolveViaTabulation");

		SolveViaTabulation();

	}

	void addPathEdge(pii e) {
		PathEdge.insert((((ull) e.fs) << 32) | e.sc);
		PathEdgeAdjRev[e.sc].pb(e.fs);
	}

	void addSummaryEdge(pii e) {
		SummaryEdge.insert((((ull) e.fs) << 32) | e.sc);
		SummaryEdgeAdj[e.fs].pb(e.sc);
	}

	void SolveViaTabulation() {   // Lines 1-12

		addPathEdge(mp(mainNodeInGExp, mainNodeInGExp));
		WorkList.insert((((ull) mainNodeInGExp) << 32) | mainNodeInGExp);

		assert(vertexTypeG[revMapGExp[mainNodeInGExp].fs] == START_VERTEX);

		// clock_t otherTime = clock();
		ForwardTabulateSLRPs();
		// timeForForwardTabulate = clock() - otherTime;
//		stopClock(Time, "ForwardTabulateSLRPs");

		// Lines 10-12, don't need to do all of this iteration if we want to answer one query
		X.assign(n_G, si());
		for (int n = 0; n < n_G; ++n) {
			int p = procOf[n];
			for (int d1 = 0; d1 <= D[p]; ++d1) {
				for (int d2 = 0; d2 <= D[p]; ++d2) {
					if (present(PathEdge, (((ull) nodeGExp[s[p]][d1]) << 32) | nodeGExp[n][d2]))
						X[n].insert(d2);
				}
			}
		}

		runTime = clock() - Time;
	}

	void Propagate(pii e) {
		assert(procOf[revMapGExp[e.fs].fs] == procOf[revMapGExp[e.sc].fs]);
		if (!present(PathEdge, (((ull) e.fs) << 32) | e.sc))
			addPathEdge(e), WorkList.insert((((ull) e.fs) << 32) | e.sc);
	}

	void ForwardTabulateSLRPs() {

		while (!WorkList.empty()) {
			ull UV = *WorkList.begin();
			WorkList.erase(WorkList.begin());
			int U = UV >> 32, V = UV & FIRST_32BIT;
			int s_p = revMapGExp[U].fs, d1 = revMapGExp[U].sc;
			int n = revMapGExp[V].fs, d2 = revMapGExp[V].sc;

			assert(vertexTypeG[s_p] == START_VERTEX);
			assert(s[procOf[s_p]] == s_p);
			assert(procOf[s_p] == procOf[n]);

			if (vertexTypeG[n] == CALL_VERTEX) {
				for (auto& W : GExp[V]) {
					if (vertexTypeG[revMapGExp[W].fs] == START_VERTEX)
						Propagate(mp(W, W));
					else {
						Propagate(mp(U, W));
					}
				}
				for (auto& W : SummaryEdgeAdj[V]) {
					Propagate(mp(U, W));
				}
			} else if (vertexTypeG[n] == EXIT_VERTEX) {
				int p = procOf[n];
				for (int i = 0; i < sz(callerEdges[p]); ++i) {
					for (auto& C : callerEdges[p][i][d1]) {
						for (auto& R : retSiteEdges[p][i][d2]) {
							if (!present(SummaryEdge, (((ull) C) << 32) | R)) {
								addSummaryEdge(mp(C, R));
								for (auto& S : PathEdgeAdjRev[C]) {
									assert(procOf[revMapGExp[S].fs] == procOf[revMapGExp[R].fs]);
									assert(vertexTypeG[revMapGExp[S].fs] == START_VERTEX);
									Propagate(mp(S, R));
								}
							}
						}
					}
				}
			} else {
				for (auto& W : GExp[V]) {
					Propagate(mp(U, W));
				}
			}
		}
	}



	void build_s() {
		s.assign(n_H, -1);

		for (int i = 0; i < n_G; ++i) {
			if (vertexTypeG[i] == START_VERTEX)
				assert(s[procOf[i]] == -1), s[procOf[i]] = i;
		}

		for (int i = 0; i < n_H; ++i) {
			assert(s[i] != -1);
		}
	}

	void buildReturnSite() {
		returnSite.assign(n_G, -1);

		for (int i = 0; i < n_G; ++i)
			if (vertexTypeG[i] == CALL_VERTEX) {
				assert(sz(G[i]) == 2);
				for (auto& v : G[i]) {
					if (vertexTypeG[v] == RETURN_SITE_VERTEX)
						assert(returnSite[i] == -1), returnSite[i] = v;
				}
				assert(returnSite[i] != -1);
				assert(procOf[i] == procOf[returnSite[i]]);
			}

	}

	bool query(int V) { // V in GExp, the V of a (U, V) query, where U is fixed to be mainNodeInGExp
		int v = revMapGExp[V].fs, d2 = revMapGExp[V].sc;
		return present(X[v], d2);
	}

};


#endif //FASTIFDS_ALGO_OLD_H
