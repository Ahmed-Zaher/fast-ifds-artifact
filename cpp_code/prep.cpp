/*
 * Input:
 *
 * Files from java_output directory, according to the format specified in "FastIFDS.java"
 *
 * Purpose:
 *
 * Do some cleaning-up on the input format so that it is in-match with the IFDS model and the PACE solvers
 * Namely, we will construct the graph representation of each edge in the supergraph according to
 * the respective analysis. The call graphs have many connected components, we'll consider procedures
 * that are in the biggest component; the remaining ones are negligibly small
 *
 *
 * Now we'll have one output for the algo of the (slightly different form):
 * {
 *      n_G m_G n_H
 *      procOf[0] procOf[1] .. procOf[n_G-1]
 *      vertexTypeG[0] vertexTypeG[1] .. vertexTypeG[n_G-1] // vertexTypeG[i] is the type of CFG node i,
 * 															// it can be a start, exit, call, or return-site
 *
 * 		u_1 v_1
 * 		u_2 v_2
 * 		   .
 * 		   .
 * 		u_{m_G} v_{m_G} // m_G edges of the form "u v" for u, v in [0, n_G): edges of the supergraph (including
 * 						// the interprocedural ones)
 *      m_H
 * 		p_1 q_1
 * 		p_2 q_2
 * 		   .
 * 		   .
 * 		p_{m_H} q_{m_G} // m_H edges of the call graph; has the form "p q" meaning p calls q (both are procedures in [0, n_H))
 * 		D[0] D[1] .. D[n_H-1]
 * 		[flow function of edge 0 in ICFG]
 * 		[flow function of edge 1 in ICFG]
 * 		   .
 * 		   .
 * 		[flow function of edge m_G-1 in ICFG]
 *  {
 *		We think of flow function of edge (u, v) as a bipartite graph (Up, Down) where:
 *			- |Up| = D[procOf[u]] + 1, nodes of Up are labeled 0..D[procOf[u]]
 *			- |Down| = D[procOf[v]] + 1, nodes of Down are labeled 0..D[procOf[v]]
 *			- this should directly correspond to the "graph representation" of the actual flow function
 *				as in the paper. 0 on either sides correspond to \textbf{0} in the papers.
 *		[flow function of edge i = (u, v) in G] has the following format:
 *		D[procOf[u]]+1 lines, j'th of them (for j in [0, D[procOf[u]]]) has the form |N(j)| N(j)[1] .. N(j)[|N(j)|]:
 *			where N(j) = the nodes in Down (which take values in [0, D[procOf[v]]]) that are connected to j in Up.
 * }
 *
 *
 * We'll also output a the appropriate formats to PACE solvers
 */

#include <bits/stdc++.h>

using namespace std;

const int OO = 1e9;
const double EPS = 1e-9;

#define ndl cout << '\n'
#define sz(v) int(v.size())
#define pb push_back
#define mp make_pair
#define fs first
#define sc second
#define present(a, x) (a.find(x) != a.end())
#ifdef LOCAL
#define db(...) ({cout << "> Line " << __LINE__  \
        << ": "; _db(#__VA_ARGS__, __VA_ARGS__);})
#else
#define db(...) true
#endif

template<class T>
void _db(const char *dbStr, T e) {
	cout << dbStr << " = " << e << endl;
}

template<class T, class... L>
void _db(const char *dbStr, T e, L... r) {
	while (*dbStr != ',') cout << *dbStr++;
	cout << " = " << e << ',';
	_db(dbStr + 1, r...);
}

template<class S, class T>
ostream &operator<<(ostream &o, const map<S, T> &v) {
	o << "[";
	int i = 0;
	for (const pair<S, T> &pr: v)
		o << (!i++ ? "" : ", ") << "{"
		  << pr.fs << " : " << pr.sc << "}";
	return o << "]";
}

template<template<class, class...> class S, class T, class... L>
ostream &operator<<(ostream &o, const S<T, L...> &v) {
	o << "[";
	int i = 0;
	for (const auto &e: v) o << (!i++ ? "" : ", ") << e;
	return o << "]";
}

template<class S, class T>
ostream &operator<<(ostream &o, const pair<S, T> &pr) {
	return o << "(" << pr.fs << ", " << pr.sc << ")";
}

ostream &operator<<(ostream &o, const string &s) {
	for (const char &c: s) o << c;
	return o;
}

template<class T> using V = vector<T>;
template<class T> using VV = V<V<T>>;
template<class T> using VVV = VV<V<T>>;
using ll = long long;
using pii = pair<int, int>;
using vi = V<int>;
using vii = V<pii>;
using vvi = VV<int>;
using mii = map<int, int>;
using umii = unordered_map<int, int>;
using si = set<int>;
using usi = unordered_set<int>;
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());



const int START_VERTEX = 0;
const int EXIT_VERTEX = 1;
const int CALL_VERTEX = 2;
const int RETURN_SITE_VERTEX = 3;


int main(int arc, char *argv[]) {
#ifdef LOCAL
	auto stTime = clock();
#endif
	ios::sync_with_stdio(false);
	cout.precision(10);
	cin.tie(0);

	V<string> fileNames;


      fileNames.pb("avrora.txt");
	fileNames.pb("sunflow.txt");
	fileNames.pb("fop.txt");
	fileNames.pb("antlr.txt");
	fileNames.pb("bloat.txt");
	fileNames.pb("chart.txt");
	fileNames.pb("eclipse.txt");
	fileNames.pb("jython.txt");
	fileNames.pb("luindex.txt");
	fileNames.pb("lusearch.txt");
	fileNames.pb("pmd.txt");
	fileNames.pb("xalan.txt");
	fileNames.pb("hsqldb.txt");



	string path = "../";

	V<string> analyses;
	analyses.pb("reachability");
	analyses.pb("uninit_var");
	analyses.pb("null_ptr");

	for (auto file: fileNames) {


		for (auto& analysis : analyses) {
			/*
			 * Just taking input
			 */
			cout << "-------------------------------------------------------------------------------------------" << endl;
			cout << file << ", " << analysis << endl;
			ifstream in(path + "java_output/" + analysis + "/" + file);
			int n_G, m_G, n_H;

			in >> n_G >> m_G >> n_H;

			db(n_G, m_G, n_H);

			vi procOf(n_G);


			for (int i = 0; i < n_G; ++i)
				in >> procOf[i], assert(procOf[i] >= 0 && procOf[i] < n_H);

			vii edgeListG;

			for (int i = 0; i < m_G; ++i) {
				int u, v;
				in >> u >> v;
				assert(u >= 0 && u < n_G);
				assert(v >= 0 && v < n_G);
				assert(u != v);
				edgeListG.pb(mp(u, v));
			}

			int m_H;

			in >> m_H;

			vii edgeListH;


			for (int i = 0; i < m_H; ++i) {
				int u, p;
				in >> u >> p;
				assert(u >= 0 && u < n_G);
				assert(p >= 0 && p < n_H);

				edgeListH.pb(mp(u, p));
			}

			m_H = sz(edgeListH);

			vi D(n_H);

			for (int i = 0; i < n_H; ++i)
				in >> D[i];

			vvi params(n_H);
			V<vii> calls(n_G);
			vi rets(n_G, -1);
			vvi assignments(n_G);

			if (analysis != "reachability") {
				for (int p = 0; p < n_H; ++p) {
					int paramsSz;
					string dummy;
					in >> dummy;
					assert(dummy == "p");
					in >> paramsSz;
					for (int i = 0; i < paramsSz; ++i) {
						int param;
						in >> param;
						assert(param >= 1 && param <= D[p]);
						params[p].pb(param);
					}
				}

				for (int u = 0; u < n_G; ++u) {
					int callsSz;
					string dummy;
					in >> dummy;
					assert(dummy == "c");
					in >> callsSz;
					for (int i = 0; i < callsSz; ++i) {
						int x, y;
						in >> x >> y;
						calls[u].pb(mp(x, y));
					}
				}

				for (int u = 0; u < n_G; ++u) {
					in >> rets[u];
					assert((rets[u] == -1) || (rets[u] >= 1 && rets[u] <= D[procOf[u]]));
				}

				for (int u = 0; u < n_G; ++u) {
					int assignmentSz;
					string dummy;
					in >> dummy;
					assert(dummy == "a");
					in >> assignmentSz;
					for (int i = 0; i < assignmentSz; ++i) {
						int x;
						in >> x;
						assert(x >= 1 && x <= D[procOf[u]]);
						assignments[u].pb(x);
					}
				}
			}

			string done;
			in >> done;

			assert(done == "done");

			/*
			 * The O represents output, and they have slightly different meaning:
			 *      - OedgeListH containts edges (p1, p2) for p1, p2 in [0, n_H) where p1 != p2
			 *          to denote that there was a call from p1 to p2. It is directed.
			 *      - OedgeListG contains all the edges of the ICFG (including inter-procedural edges)
			 *          has edges of the form (u, v) for u, v in [0, n_G) to denote flow from u to v.
			 *      - The rest of the stuff have the same meaning
			 */

			V<vvi> flow;
			vvi adjCG(n_H), adjCGRev(n_H);
			vvi adjCGUndir(n_H);
			vi vertexTypeG;

			{
				int On_G = n_G, Om_G = 0;
				vii OedgeListG;
				set<pii> OedgeListH;


				vi inDeg(n_G);
				vi outDeg(n_G);


				for (auto pr : edgeListG) {
					int u = pr.fs, v = pr.sc;
					assert(procOf[u] == procOf[v]);
					++outDeg[u], ++inDeg[v];
				}

				vvi starts(n_H), ends(n_H);


				for (int u = 0; u < On_G; u++) {
					assert(procOf[u] >= 0 && procOf[u] < n_H);
					if (!inDeg[u])
						starts[procOf[u]].pb(u);
					if (!outDeg[u])
						ends[procOf[u]].pb(u);
				}


				vi superStart(n_H), superEnd(n_H);
				vector<int> callRet(On_G, -1);

				// IDs of superstarts and superends start after n_G
				for (int p = 0; p < n_H; ++p) {

					assert(sz(starts[p]) && sz(ends[p]));

					superStart[p] = On_G++;
					procOf.pb(p);
					callRet.pb(-1);

					superEnd[p] = On_G++;
					procOf.pb(p);
					callRet.pb(-1);

					for (int st : starts[p]) {
						// Case 1. super-start to start: \S. S union (Loc[p] \ Params[p])
						int e = sz(OedgeListG);
						++Om_G;
						OedgeListG.pb(mp(superStart[p], st));
						flow.pb(vvi(D[p] + 1, vi()));


						if (analysis == "reachability") {
							flow[e][0].pb(0);
							flow[e][1].pb(1);
						} else {
							si complementProcParams;

							for (int lcl = 1; lcl <= D[p]; ++lcl)
								complementProcParams.insert(lcl);

							for (auto& elem : params[p])
								complementProcParams.erase(elem);

							flow[e][0].pb(0);

							for (auto& elem : complementProcParams)
								flow[e][0].pb(elem);

							for (int lcl = 1; lcl <= D[p]; ++lcl) if (!present(complementProcParams, lcl))
									flow[e][lcl].pb(lcl);
						}
					}


					for (int en : ends[p]) {
						// Case 2. end to super-end: \S. S
						int e = sz(OedgeListG);
						++Om_G;
						OedgeListG.pb(mp(en, superEnd[p]));

						flow.pb(vvi(D[p] + 1, vi()));

						for (int lcl = 0; lcl <= D[p]; ++lcl)
							flow[e][lcl].pb(lcl);
					}
				}

				// ID for invoke = ID for call node, the ID for return-site will takes values after n_G+#superstarts+#superends
				// callRet[i] = -1 if i is not invoke node, otherwise it's the return-site nodes
				for (auto edge : edgeListH) {
					int u = edge.fs;
					assert(callRet[u] == -1);
					callRet[u] = On_G++;
					procOf.pb(procOf[u]);
					callRet.pb(-1);
				}

				for (auto edge : edgeListH) {
					int u = edge.fs, p = edge.sc;

					// Case 3. call-to-return-site:
					// we assume call-by-ref:   \S. S \ Params[procOf[u]]

					int e = sz(OedgeListG);
					++Om_G;
					OedgeListG.pb(mp(u, callRet[u]));
					flow.pb(vvi(D[procOf[u]] + 1, vi()));
					flow[e][0].pb(0);


					if (analysis != "reachability") {


						si complementProcParams;

						for (int lcl = 1; lcl <= D[procOf[u]]; ++lcl)
							complementProcParams.insert(lcl);

						for (auto& elem : params[procOf[u]])
							complementProcParams.erase(elem);

						for (auto& elem : complementProcParams)
							flow[e][elem].pb(elem);

					}


					// Case 4. call-to-start: \S. (S intersect args)[map labels to match p]
					e = sz(OedgeListG);
					++Om_G;
					OedgeListG.pb(mp(u, superStart[p]));
					flow.pb(vvi(D[procOf[u]] + 1, vi()));

					flow[e][0].pb(0);
					if (analysis != "reachability") {
						for (auto& pr : calls[u]) {
							int lcl_1 = pr.fs, lcl_2 = pr.sc;
							assert(lcl_1 >= 1 && lcl_1 <= D[procOf[u]]);
							assert(lcl_2 >= 1 && lcl_2 <= D[p]);
							flow[e][lcl_1].pb(lcl_2);
						}
					} else {
						flow[e][1].pb(1);
					}


					// Case 5. exit-to-return-site
					// call-by-ref: \S. (S intersect params)[map labels to match caller-proc]
					e = sz(OedgeListG);
					++Om_G;
					OedgeListG.pb(mp(superEnd[p], callRet[u]));
					flow.pb(vvi(D[p] + 1, vi()));

					flow[e][0].pb(0);

					if (analysis != "reachability") {

						for (auto& pr : calls[u]) {
							int lcl_1 = pr.fs, lcl_2 = pr.sc;
							assert(lcl_1 >= 1 && lcl_1 <= D[procOf[u]]);
							assert(lcl_2 >= 1 && lcl_2 <= D[p]);
							flow[e][lcl_2].pb(lcl_1);
						}
					} else {
						flow[e][1].pb(1);
					}
				}


				for (auto edge : edgeListG) {
					int u = edge.fs, v = edge.sc;

					int e = sz(OedgeListG);
					++Om_G;
					flow.pb(vvi(D[procOf[u]] + 1, vi()));
					flow[e][0].pb(0);

					// Case 6. we have an assignment x := ...
					// then we have:
					//      \S. if (S intersect Loc[RHS]) =/= phi then: S union {x} else: S \ {x}
					if (!assignments[u].empty()) {
						int x = assignments[u][0];
						assert(x >= 1 && x <= D[procOf[u]]);

						si LocRHS;

						for (int i = 1; i < sz(assignments[u]); ++i) {
							LocRHS.insert(assignments[u][i]);
						}

						for (int i = 1; i <= D[procOf[u]]; ++i) {
							if (present(LocRHS, i)) {
								flow[e][i].pb(x);
								if (i != x)
									flow[e][i].pb(i);
							} else {
								if (i != x)
									flow[e][i].pb(i);
							}
						}
					} else {
						// Case 7. None of the above cases, just use \S. S, reachability will always get to this branch
						for (int i = 1; i <= D[procOf[u]]; ++i)
							flow[e][i].pb(i);
					}

					if (callRet[u] != -1)
						OedgeListG.pb(mp(callRet[u], v));
					else
						OedgeListG.pb(mp(u, v));
				}

				assert(sz(flow) == sz(OedgeListG));



				for (auto edge : edgeListH) {
					int u = edge.fs, p = edge.sc;
					if (procOf[u] != p) {
						adjCG[procOf[u]].pb(p);
						adjCGRev[p].pb(procOf[u]);
						OedgeListH.insert(mp(procOf[u], p));

						adjCGUndir[procOf[u]].pb(p);
						adjCGUndir[p].pb(procOf[u]);
					}
				}


				// generate vertexTypeG

				assert(On_G == sz(procOf));
				assert(Om_G == sz(OedgeListG));

				vertexTypeG.assign(On_G, -1);

				for (int i = 0; i < n_H; ++i)
					vertexTypeG[superStart[i]] = START_VERTEX, vertexTypeG[superEnd[i]] = EXIT_VERTEX;

				for (int i = 0; i < On_G; ++i)
					if (callRet[i] != -1)
						vertexTypeG[i] = CALL_VERTEX, vertexTypeG[callRet[i]] = RETURN_SITE_VERTEX;

				// Now edgeListG, edgeListH, n_G, m_G are updated
				n_G = On_G, m_G = Om_G;
				edgeListG = OedgeListG;

				edgeListH.clear();

				for (auto& pr : OedgeListH)
					edgeListH.pb(pr);


			}

			// Transforming the problem to include only biggest component:
			{

				int biggestCompClr = -1;
				vi methodCCClr(n_H);    // methodCCClr[p] = methodCCClr[p'] if p and p' are in the same CC in the
				si keep; // set of procedures that we'll keep
				// undirected version of the CG
				{
					int ccClr = 0, ccSz = 0;
					vii ccSizes; // stores pairs (size of CC, ccID)
					mii ccSizesFreq;
					vi vis(n_H);
					function<void(int)> dfs = [&](int u) {
						vis[u] = true;
						++ccSz;
						methodCCClr[u] = ccClr;
						for (auto& v : adjCGUndir[u]) {
							if (!vis[v])
								dfs(v);
						}
					};
					for (int i = 0; i < n_H; ++i) {
						if (!vis[i]) {
							ccSz = 0;
							dfs(i);
							ccSizes.pb(mp(ccSz, ccClr));
							ccSizesFreq[ccSz]++;
							ccClr++;
						}
					}
					db(ccSizesFreq);
					sort(ccSizes.rbegin(), ccSizes.rend());
					biggestCompClr = ccSizes[0].sc;

					for (int i = 0; i < n_H; ++i) {
						if (methodCCClr[i] == biggestCompClr)
							keep.insert(i);
					}
					int  biggestCompSz = ccSizes[0].fs;

				}


				assert(biggestCompClr != -1);


				// OO stands for instance considering only the trimmed call graph, after we build it we set data := Tdata
				int Tn_G, Tn_H;
				vi TprocOf;
				vii TedgeListG;
				vii TedgeListH;
				vi TvertexTypeG;
				vi TD;
				V<vvi> Tflow;

				vi ICFGMapping(n_G, -1);   // relabeling ICFG nodes to include only those in methods in biggest comp
				vi methodMapping(n_H, -1); //  relabeling CG nodes to include only those methods in biggest comp

				Tn_G = Tn_H = 0;

				for (int i = 0; i < n_H; ++i) {
					if (present(keep, i))
						methodMapping[i] = Tn_H++;
				}

				for (int i = 0; i < n_G; ++i) {
					if (present(keep, procOf[i]))
						ICFGMapping[i] = Tn_G++;
				}


				TprocOf.assign(Tn_G, -1);

				for (int i = 0; i < n_G; ++i) if (ICFGMapping[i] != -1) {
					TprocOf[ICFGMapping[i]] = methodMapping[procOf[i]];
				}

				for (int i = 0; i < Tn_G; ++i) {
					assert(TprocOf[i] != -1);
					assert(TprocOf[i] < Tn_H);
				}

				int eId = 0;

				for (auto& pr : edgeListG) {
					int u = pr.fs, v = pr.sc;
					int sat = (ICFGMapping[u] != -1) + (ICFGMapping[v] != -1);
					assert(sat == 0 || sat == 2);   // holds if we take whole comp
					if (sat == 2)
						TedgeListG.pb(mp(ICFGMapping[u], ICFGMapping[v])), Tflow.pb(flow[eId]);
					++eId;
				}

				for (auto& pr : edgeListH) {
					int p1 = pr.fs, p2 = pr.sc;
					int sat = (methodMapping[p1] != -1) + (methodMapping[p2] != -1);
					assert(sat == 0 || sat == 2);   // holds if we take whole comp
					if (sat == 2)
						TedgeListH.pb(mp(methodMapping[p1], methodMapping[p2]));
				}


				TvertexTypeG.assign(Tn_G, -2);

				for (int i = 0; i < n_G; ++i) if (ICFGMapping[i] != -1) {
					TvertexTypeG[ICFGMapping[i]] = vertexTypeG[i];
				}



				for (int i = 0; i < Tn_G; ++i) {
					assert(TvertexTypeG[i] != -2);
				}

				TD.assign(Tn_H, -1);

				for (int i = 0; i < n_H; ++i) if (methodMapping[i] != -1) {
						TD[methodMapping[i]] = D[i];
					}

				for (int i = 0; i < Tn_H; ++i) {
					assert(TD[i] != -1);
				}

				n_G = Tn_G, n_H = Tn_H;
				procOf = TprocOf;
				edgeListG = TedgeListG;
				edgeListH = TedgeListH;
				vertexTypeG = TvertexTypeG;
				D = TD;
				flow = Tflow;
			}


			// Printing output
			{
				// 1-based input to be fed to the treedepth PACE solver, make it 0-based in IfdsIntance
				// same over all analyses, will just overwrite the same thing, no problem

				set<pii> PACEedges;
				for (auto edge : edgeListH) {
					int u = edge.fs, v = edge.sc;
					PACEedges.insert(mp(min(u, v), max(u, v)));
				}

				ofstream out(path + "treedepth_solver_input/" + file);
				out << "p tdp " << n_H << ' ' << sz(PACEedges) << '\n';

				// 1-based connected graph with no self loops of multiple edges to fit the PACE solver

				for (auto pr : PACEedges) {
					out << pr.fs + 1 << ' ' << pr.sc + 1 << '\n';
				}

				{
					vvi adjPACE(n_H + 1);
					for (auto& edge : PACEedges) {
						int u = edge.fs + 1, v = edge.sc + 1;
						if (u == v)
							db(u, v);
						assert(u != v);
						assert(u >= 1 && u <= n_H);
						assert(v >= 1 && v <= n_H);
						adjPACE[u].pb(v);
						adjPACE[v].pb(u);
					}
					int cnt = 0;
					vi vis(n_H + 1);
					function<void(int)> dfs = [&](int u) {
						vis[u] = true;
						++cnt;
						for (auto& v : adjPACE[u]) {
							if (!vis[v])
								dfs(v);
						}
					};
					dfs(1);
					assert(cnt == n_H);
				}
			}

			;
			// IMP: uncomment
			{
				// 1-based input to be fed to the treewidth PACE solver, make it 0-based in IfdsIntance
				// same over all analyses, will just overwrite the same thing, no problem

				ofstream out(path + "treewidth_solver_input/" + file);

				set<pii> treewidthEdges;

				for (auto& edge : edgeListG) if (procOf[edge.fs] == procOf[edge.sc])
					treewidthEdges.insert(mp(min(edge.fs, edge.sc), max(edge.fs, edge.sc)));

				int diff = 0;
				for (int i = 0; i < n_G; ++i) if (vertexTypeG[i] == START_VERTEX)
						treewidthEdges.insert(mp(i, n_G)), ++diff;

				assert(diff == n_H);

				out << "p tw " << n_G + 1 << ' ' << sz(treewidthEdges) << '\n';

				// 1-based connected graph with no self loops of multiple edges to fit the PACE solver

				for (auto pr : treewidthEdges) {
					out << pr.fs + 1 << ' ' << pr.sc + 1 << '\n';
				}

				{
					vvi adjPACE(n_G + 2);
					for (auto& edge : treewidthEdges) {
						int u = edge.fs + 1, v = edge.sc + 1;
						if (u == v)
							db(u, v);
						assert(u != v);
						assert(u >= 1 && u <= n_G + 1);
						assert(v >= 1 && v <= n_G + 1);
						adjPACE[u].pb(v);
						adjPACE[v].pb(u);
					}
					int cnt = 0;
					vi vis(n_G + 2);
					function<void(int)> dfs = [&](int u) {
						vis[u] = true;
						++cnt;
						for (auto& v : adjPACE[u]) {
							if (!vis[v])
								dfs(v);
						}
					};
					dfs(1);
					assert(cnt == n_G + 1);
				}
			}


			ofstream out(path + "prep_output/" + analysis + "/" + file);

			// Print results according to input format of main.cpp

			out << n_G << ' ' << sz(edgeListG) << ' ' << n_H << '\n';

			for (int i = 0; i < n_G; ++i)
				out << procOf[i] << " \n"[i + 1 == n_G];

			for (int i = 0; i < n_G; ++i)
				out << vertexTypeG[i] << " \n"[i + 1 == n_G];

			for (auto pr : edgeListG)
				out << pr.fs << ' ' << pr.sc << '\n';

			out << sz(edgeListH) << '\n';

			for (auto pr : edgeListH)
				out << pr.fs << ' ' << pr.sc << '\n';


			for (int p = 0; p < n_H; ++p)
				out << D[p] << " \n"[p + 1 == n_H];


			assert(sz(edgeListG) == sz(flow));

			m_G = sz(edgeListG);

			for (int i = 0; i < m_G; ++i) {
				int u = edgeListG[i].fs, v = edgeListG[i].sc;
				assert(sz(flow[i]) == D[procOf[u]] + 1);
				for (int j = 0; j <= D[procOf[u]]; ++j) {
					out << sz(flow[i][j]);
					for (auto& k : flow[i][j]) {
					    assert(k >= 0 && k <= D[procOf[v]]);
						out << " " << k;
					}
					out << "\n";
				}
			}

			out << "done\n";

			db(n_G, m_G, n_H);

		}


	}

	return 0;
}
