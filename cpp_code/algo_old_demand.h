#ifndef FASTIFDS_ALGO_OLD_DEMAND_H
#define FASTIFDS_ALGO_OLD_DEMAND_H

#include "template.h"
#include "SimpleIfdsInstance.h"

int HashModG_alg_old_demand;



class algo_old_demand {
public:
	void stopClock(double Time, string s = "") {

		cout << "[" << s << "]" << " time taken = " <<
		     (omp_get_wtime() - Time) << " s" << endl;
	}
	/*
	 * We'll do the implementation according to the conference version
	 * The U is assumed to be of the form (s_p, 0) for some p in [0, n_H).
	 */
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

	int mainNodeInGExp; // The U in a (U, V) queryUisStart, has the form (s_p, 0)
	vi s;               // s[p] for p in [0, n_H) is a node in G which is the start of proc p.
	vi callerNode;      // callerNode[u] for u in [0, n_G) and u is a RETURN_SITE_VERTEX, is the corresponding
						// CALL_VERTEX in G.

	// The rest is direct translation of the paper, line numbering is according to the original version
	// (we make the variables global instead of passing them all the time as in the conference-version)

	unordered_set<ull> PathEdge, SummaryEdge;
	unordered_set<int> ReachableNodes, VisitedNodes, ReachableNodesRelevantToDemand;
	unordered_set<int> ReachableNodesFromSrc;   // auxiliary
	vvi SummaryEdgeAdj, SummaryEdgeAdjRev, PathEdgeAdj;
	vvi GExpRev;
	V<V<V<V<int>>>> callerEdges, retSiteEdges;

	ull FIRST_32BIT;

	double timeForIsMemberOfSolution;
	double runTime;
	double Time;





	algo_old_demand(SimpleIfdsInstance* instance, int _mainNodeInGExp)
			: mainNodeInGExp(_mainNodeInGExp) {

		Time = omp_get_wtime();

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

		HashModG_alg_old_demand = n_GExp;
		buildCallerNode();



		SummaryEdgeAdjRev.assign(n_GExp, vi());
		SummaryEdgeAdj.assign(n_GExp, vi());
		PathEdgeAdj.assign(n_GExp, vi());

		// Building GExpRev
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

		for (int r = 0; r < n_G; ++r) {
			if (vertexTypeG[r] == RETURN_SITE_VERTEX) {
				int c = callerNode[r];
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

		ReachableNodes.insert(mainNodeInGExp);

		runTime = omp_get_wtime() - Time;
	}

	void addPathEdge(pii e) {
		PathEdge.insert((((ull) e.fs) << 32) | e.sc);
		PathEdgeAdj[e.fs].pb(e.sc);
	}

	void addSummaryEdge(pii e) {
		SummaryEdge.insert((((ull) e.fs) << 32) | e.sc);
		SummaryEdgeAdjRev[e.sc].pb(e.fs);
		SummaryEdgeAdj[e.fs].pb(e.sc);
	}

	bool IsMemberOfSolution(int ndash, int ddash) {
		// the declare line
		unordered_set<int> NodeWorkList;
		unordered_set<ull> EdgeWorkList;

		// Line 1-2
		ReachableNodesRelevantToDemand.clear();
		Visit(nodeGExp[ndash][ddash], NodeWorkList);


		while (!NodeWorkList.empty()) {
			int N = *NodeWorkList.begin();
			NodeWorkList.erase(NodeWorkList.begin());
			int n = revMapGExp[N].fs, d = revMapGExp[N].sc;

			if (vertexTypeG[n] == RETURN_SITE_VERTEX) {
				int c = callerNode[n];
				assert(sz(G[c]) == 2);
				int p = -1;
				for (auto& elem : G[c]) {
				    if (vertexTypeG[elem] == START_VERTEX)
						assert(p == -1), p = procOf[elem];
				}
				assert(p != -1);

				EdgeWorkList.clear();

				for (auto& W : GExpRev[N]) {
					if (vertexTypeG[revMapGExp[W].fs] == EXIT_VERTEX)
						Propagate(mp(W, W), EdgeWorkList);
				}
				BackwardTabulateSLRPs(EdgeWorkList);
				for (auto& C : GExpRev[N]) {
				    if (vertexTypeG[revMapGExp[C].fs] == EXIT_VERTEX)
						Visit(C, NodeWorkList);
				}
				for (auto& C : SummaryEdgeAdjRev[N]) {
					assert(vertexTypeG[revMapGExp[C].fs] == CALL_VERTEX);
					Visit(C, NodeWorkList);
				}

			} else if (vertexTypeG[n] == START_VERTEX) {
				for (int i = 0; i < sz(callerEdges[procOf[n]]); ++i) {
					for (auto& U : callerEdges[procOf[n]][i][d]) {
					    assert(vertexTypeG[revMapGExp[U].fs] == CALL_VERTEX), Visit(U, NodeWorkList);
					}
				}
			} else {
				for (auto& V : GExpRev[N]) {
				    Visit(V, NodeWorkList);
				}
			}
		}
		UpdateReachableNodes();
		return present(ReachableNodes, nodeGExp[ndash][ddash]);
	}

	void Visit(int N, unordered_set<int>& NodeWorkList) {
		if (present(ReachableNodes, N))
			ReachableNodesRelevantToDemand.insert(N);
		else if (!present(VisitedNodes, N))
			VisitedNodes.insert(N), NodeWorkList.insert(N);
	}

	void UpdateReachableNodes() {
		while (!ReachableNodesRelevantToDemand.empty()) {
			int N = *ReachableNodesRelevantToDemand.begin();
			ReachableNodesRelevantToDemand.erase(ReachableNodesRelevantToDemand.begin());

			for (auto& M : GExp[N]) {

				if ((vertexTypeG[revMapGExp[N].fs] == EXIT_VERTEX)
				 && (vertexTypeG[revMapGExp[M].fs] == RETURN_SITE_VERTEX))
					continue;
				if (present(VisitedNodes, M) && !present(ReachableNodes, M))
					ReachableNodes.insert(M), ReachableNodesRelevantToDemand.insert(M);
			}

			for (auto& M : SummaryEdgeAdj[N]) {

				if ((vertexTypeG[revMapGExp[N].fs] == EXIT_VERTEX)
				    && (vertexTypeG[revMapGExp[M].fs] == RETURN_SITE_VERTEX))
					continue;
				if (present(VisitedNodes, M) && !present(ReachableNodes, M))
					ReachableNodes.insert(M), ReachableNodesRelevantToDemand.insert(M);
			}
		}
	}

	void BackwardTabulateSLRPs(unordered_set<ull>& EdgeWorkList) {

		while (!EdgeWorkList.empty()) {
			ull UV = *EdgeWorkList.begin();
			EdgeWorkList.erase(EdgeWorkList.begin());
			int U = UV >> 32, V = UV & FIRST_32BIT;
			int n = revMapGExp[U].fs, d2 = revMapGExp[U].sc;
			int e_p = revMapGExp[V].fs, d1 = revMapGExp[V].sc;

			assert(vertexTypeG[e_p] == EXIT_VERTEX);

			if (vertexTypeG[n] == RETURN_SITE_VERTEX) {
				for (auto& W : GExpRev[U]) {
				    if (vertexTypeG[revMapGExp[W].fs] == EXIT_VERTEX)
					    Propagate(mp(W, W), EdgeWorkList);
					else {
						if (W != mainNodeInGExp)
							assert(vertexTypeG[revMapGExp[W].fs] == CALL_VERTEX);
					    Propagate(mp(W, V), EdgeWorkList);
					}
				}
				for (auto& C : SummaryEdgeAdjRev[U]) {
					Propagate(mp(C, V), EdgeWorkList);
				}
			} else if (vertexTypeG[n] == START_VERTEX) {
				int p = procOf[n];
				for (int i = 0; i < sz(callerEdges[p]); ++i) {
					for (auto& C : callerEdges[p][i][d2]) {
					    for (auto& R : retSiteEdges[p][i][d1]) {
					        int c = revMapGExp[C].fs, r = revMapGExp[R].fs;
							assert(procOf[c] == procOf[r]);
							assert(callerNode[r] == c);
							if (!present(SummaryEdge, (((ull) C) << 32) | R)) {
								addSummaryEdge(mp(C, R));
								for (auto& E : PathEdgeAdj[R]) {
									Propagate(mp(C, E), EdgeWorkList);
								}
							}
					    }
					}
				}
			} else {
				for (auto& W : GExpRev[U]) {
					Propagate(mp(W, V), EdgeWorkList);
				}
			}
		}
	}

	void Propagate(pii e, unordered_set<ull>& EdgeWorkList) {
		if (!present(PathEdge, (((ull) e.fs) << 32) | e.sc))
			addPathEdge(e), EdgeWorkList.insert((((ull) e.fs) << 32) | e.sc);
	}




	void buildCallerNode() {
		callerNode.assign(n_G, -1);

		for (int i = 0; i < n_G; ++i)
			if (vertexTypeG[i] == CALL_VERTEX) {
				assert(sz(G[i]) == 2);
				for (auto& v : G[i])
					if (vertexTypeG[v] == RETURN_SITE_VERTEX)
						assert(callerNode[v] == -1), callerNode[v] = i;
			}

		for (int i = 0; i < n_G; ++i) {
			if (vertexTypeG[i] == RETURN_SITE_VERTEX)
				assert(callerNode[i] != -1), assert(procOf[i] == procOf[callerNode[i]]);
		}

	}


	bool query(int V) { // V in GExp, the V of a (U, V) queryUisStart, where U is fixed to be mainNodeInGExp

		Time = omp_get_wtime();
		bool ans = IsMemberOfSolution(revMapGExp[V].fs, revMapGExp[V].sc);
		timeForIsMemberOfSolution = omp_get_wtime() - Time;
		return ans;
	}

};



#endif //FASTIFDS_ALGO_OLD_DEMAND_H
