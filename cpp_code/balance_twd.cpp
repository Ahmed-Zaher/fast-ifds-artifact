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

/*
 * Input:
 *
 * Unbaalanced treewidth decomp from the PACE solver in treewidth_solver_output and some relevant info from prep_output
 *
 * Purpose:
 *
 * Balance the treewidth decompositions and store them in balanced_twd
 *
 * Next phases:
 *
 * The FastIFDS algorithm
 *
 * Output:
 *
 * Will output treewidth-decomposition of the following format
 *      twd 0
 *      [treewidth decomp. of proc 0]
 *      twd 1
 *      [treewidth decomp. of proc 1]
 *          .
 *          .
 *      twd n_H - 1
 *      [treewidth decomp. of proc n_H - 1]
 *   {
 *		[treewidth decomp. of proc i] has the following format:
 *
 *		n_B (number of bags, the bags are labeled with [0, n_B))
 *		n_B lines of the form:
 *			bagSz [bagSz numbers], the elements in the bag with their indices in G
 *		n_B space-separated lines, i'th num (i is 0-based) is the parent of bag i - root has parent = -1
 *	}
 */




const int START_VERTEX = 0;
const int EXIT_VERTEX = 1;
const int CALL_VERTEX = 2;
const int RETURN_SITE_VERTEX = 3;

void validateTWD(VVV<int> TWD, vi TWD_root, vvi TWD_par, vi procOf, vii edgeListG) {
	int n_G = sz(procOf), n_H = sz(TWD);
	vector<vector<int>> bagsContainingNode(n_G);

	VVV<int> TWD_adj(n_H);

	for (int p = 0; p < n_H; ++p) {

		assert(sz(TWD[p]) == sz(TWD_par[p]));
		int n_T_i = sz(TWD[p]); // number of bags in treew. dec. of G_i
		TWD_adj[p].assign(n_T_i, vi());

		for (int j = 0; j < n_T_i; ++j)
			if (TWD_par[p][j] != -1) {
				TWD_adj[p][j].pb(TWD_par[p][j]);
				TWD_adj[p][TWD_par[p][j]].pb(j);
			}
	}

	vvi twd_Depth(n_H);
	vi rb(n_G, -1);

	int maxDepth = 0, biggestBagSz = 0;
	double avgBiggestBagSz = 0;

	function<void(vvi&, int, int, int, int)> dfsRLocal = [&] (vvi& adj, int u, int parent, int p, int dep) {

		maxDepth = max(maxDepth, dep);

		for (auto& v : TWD[p][u]) {
			if (rb[v] == -1)
				rb[v] = u;
			bagsContainingNode[v].pb(u);
		}

		twd_Depth[p][u] = dep;
		for (int v : adj[u])
			if (v != parent)
				dfsRLocal(adj, v, u, p, dep + 1);
	};


	for (int p = 0; p < n_H; ++p) {
		twd_Depth[p].assign(sz(TWD[p]), -1);
		dfsRLocal(TWD_adj[p], TWD_root[p], -1, p, 0);
		for (auto& elem : twd_Depth[p]) {
			assert(elem >= 0);
		}
		int procBiggestBagSz = 0;
		for (auto& elem : TWD[p]) {
			biggestBagSz = max(biggestBagSz, sz(elem));
			procBiggestBagSz = max(procBiggestBagSz, sz(elem));
		}
		avgBiggestBagSz += procBiggestBagSz;
	}

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
//		db(u);
		for (auto& b : bagsContainingNode[u]) {
//		    db(b, twd_Depth[p][b]);
		}
		for (auto& b : bagsContainingNode[u]) if (b != rb[u]) {
//			db(TWD_par[p][b], twd_Depth[p][b], twd_Depth[p][TWD_par[p][b]]);
				assert(present(prevBags, TWD_par[p][b]));
				prevBags.insert(b);
			}
	}

	// check that each edge appears in a bag

	vvi G(n_G);

	for (auto& edge : edgeListG) {
		G[edge.fs].pb(edge.sc);
	}

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

	cout << "Valid TWD!!" << endl;
	db(maxDepth, biggestBagSz, avgBiggestBagSz / n_H);
}

void fixTWD(VVV<int>& TWD, vi TWD_root, vvi TWD_par, vi procOf, vii edgeListG) {
	int n_G = sz(procOf), n_H = sz(TWD);
	vector<vector<int>> bagsContainingNode(n_G);

	VVV<int> TWD_adj(n_H);

	for (int p = 0; p < n_H; ++p) {

		assert(sz(TWD[p]) == sz(TWD_par[p]));
		int n_T_i = sz(TWD[p]); // number of bags in treew. dec. of G_i
		TWD_adj[p].assign(n_T_i, vi());

		for (int j = 0; j < n_T_i; ++j)
			if (TWD_par[p][j] != -1) {
				TWD_adj[p][j].pb(TWD_par[p][j]);
				TWD_adj[p][TWD_par[p][j]].pb(j);
			}
	}

	vvi twd_Depth(n_H);
	vi rb(n_G, -1);

	int cnt;

	function<void(vvi&, int, int, int, int)> dfsRLocal = [&] (vvi& adj, int u, int parent, int p, int dep) {
		++cnt;
		for (auto& v : TWD[p][u]) {
			if (rb[v] == -1)
				rb[v] = u;
			bagsContainingNode[v].pb(u);
		}


		twd_Depth[p][u] = dep;
		for (int v : adj[u])
			if (v != parent)
				dfsRLocal(adj, v, u, p, dep + 1);
	};

	vi mxBagSz(n_H, 0);


	for (int p = 0; p < n_H; ++p) {
		twd_Depth[p].assign(sz(TWD[p]), -1);

		cnt = 0;
		dfsRLocal(TWD_adj[p], TWD_root[p], -1, p, 0);
//		assert(cnt == sz(TWD[p]));   // make sure TWD is connected

		for (auto& bag : TWD[p]) {
			mxBagSz[p] = max(mxBagSz[p], sz(bag));
		}
	}


	for (int u = 0; u < n_G; ++u) {
		sort(bagsContainingNode[u].begin(), bagsContainingNode[u].end(), [&](int b1, int b2) {
			return twd_Depth[procOf[u]][b1] < twd_Depth[procOf[u]][b2];
		});
	}



	// check that bagsContainingNode[u] is connected and non-empty

	int mxDiff = 0;
	int totDiff = 0;


	for (int u = 0; u < n_G; ++u) {
		si prevBags;
		assert(!bagsContainingNode[u].empty());
		assert(rb[u] != -1);
		assert(bagsContainingNode[u][0] == rb[u]);
		prevBags.insert(rb[u]);
		int p = procOf[u];

		for (auto& b : bagsContainingNode[u]) if (b != rb[u]) {
//			db(TWD_par[p][b], twd_Depth[p][b], twd_Depth[p][TWD_par[p][b]]);'
				assert(twd_Depth[p][b]);
				int curBag = TWD_par[p][b];
				while (!present(prevBags, curBag))
					prevBags.insert(curBag), TWD[p][curBag].pb(u), curBag = TWD_par[p][curBag];
				prevBags.insert(b);
			}
		mxDiff = max(mxDiff, sz(prevBags) - sz(bagsContainingNode[u]));
		totDiff += sz(prevBags) - sz(bagsContainingNode[u]);
	}

	vi mxBagSzNew(n_H, 0);

	int mxDiffTw = 0, totDiffTw = 0;
	double mxDiffToOrig = 0;

	for (int p = 0; p < n_H; ++p) {
		for (auto& bag : TWD[p]) {
			mxBagSzNew[p] = max(mxBagSzNew[p], sz(bag));
		}
		mxDiffTw = max(mxDiffTw, mxBagSzNew[p] - mxBagSz[p]);
		mxDiffToOrig = max(mxDiffToOrig, double(mxBagSzNew[p] - mxBagSz[p]) / mxBagSz[p]);
		totDiffTw += mxBagSzNew[p] - mxBagSz[p];
	}



	db(mxDiff, totDiff, double(totDiff) / n_G);
	db(mxDiffTw, totDiffTw, double(totDiffTw) / n_H);
	db(mxDiffToOrig);

	// check that each edge appears in a bag

	vvi G(n_G);

	for (auto& edge : edgeListG) {
		G[edge.fs].pb(edge.sc);
	}

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

	cout << "fixed TWD!!" << endl;
}

void rmv(V<si>& adj, int u, int v) {
	assert(min(u, v) >= 0);
	assert(max(u, v) < sz(adj));
	adj[u].erase(v);
	adj[v].erase(u);
}


void add(V<si>& adj, int u, int v) {
	assert(min(u, v) >= 0);
	assert(max(u, v) < sz(adj));
	adj[u].insert(v);
	adj[v].insert(u);
}

pair<pair<V<si>, int>, vvi> getTWDofT(V<si> adj_T0, int root_T0) {
	// Given a binary tree adj_T, find a treewidth-decomp.
	// of the TREE ITSELF that has O(log(n_T)) depth and width 2
	// ret.fs.fs = the tree, ret.fs.sc = root, ret.sc = bags

	int n = sz(adj_T0);

	V<pair<V<si>, int>> T;

	V<si> TDash(n);
	vi par_TDash(n, -1);
	int TDashRoot = -1;
	vvi TDash_bags;
	V<si> B;
	V<si> N;

	for (int i = 0; i < n; ++i) {
		B.pb(si({i}));
		N.pb(adj_T0[i]);
		TDash_bags.pb(vi({i}));
	}

	T.pb(mp(adj_T0, root_T0));


	function<int(V<si>)> cntNonEmpty = [](V<si> x) {
		int ret = 0;
		for (auto& y : x) {
		    ret += !y.empty();
		}
		return ret;
	};

	while (cntNonEmpty(T.back().fs) > 1) {

//		db(sz(T.back().fs));

		V<si> adj = T.back().fs;
		int root = T.back().sc;


		assert(root >= 0 && root < n);

		vi par(n, -1), dep(n, 0);

		function<void(int, int, int)> dfs = [&](int u, int p, int d) {
//			db(u, p, d);
			assert(u >= 0 && u < n);
			par[u] = p;
			dep[u] = d;

			for (auto& v : adj[u]) if (v != p) {
			    dfs(v, u, d + 1);
			}
		};

		dfs(root, -1, 0);


		vi q;

		for (int i = 0; i < n; ++i) if (i != root && !adj[i].empty()) {
			q.pb(i);
		}



		sort(q.begin(), q.end(), [&](int u, int v) { return dep[u] > dep[v]; });

		V<bool> taken(n);

		V<si> adj_new = adj;
		int root_new = root;

		for (auto& u : q) {
			assert(par[u] != -1);

			int v = par[u];

			if (!taken[u] && !taken[v] && min(sz(adj[u]), sz(adj[v])) <= 2) {
				if (mp(sz(adj[v]), sz(adj[u])) != mp(3, 2)) {
					int w = n++;

					adj_new.pb(si());
					TDash.pb(si());
					par_TDash.pb(-1);
					TDash_bags.pb(vi());
					B.pb(si());
					N.pb(si());

					par_TDash[u] = par_TDash[v] = w;

					add(TDash, u, w);
					add(TDash, v, w);

					taken[u] = taken[v] = 1;


					si N_u = adj_new[u], N_v = adj_new[v];

					for (auto& x : N_u) {
						rmv(adj_new, u, x);
						if (x != v) {
						    add(adj_new, w, x);
						}
					}

					for (auto& x : N_v) {
						if (x != u) {
							rmv(adj_new, v, x);
							add(adj_new, w, x);
						}
					}

					if (v == root_new)
						root_new = par_TDash[root_new];



					for (auto& x : B[u]) {
						bool include = false;
						for (auto& v : adj_T0[x]) {

						}
					}

					for (auto& x : N[u]) if (!present(B[v], x)) N[w].insert(x);
					for (auto& x : N[v]) if (!present(B[u], x)) N[w].insert(x);

					for (auto& x : B[u]) {
						for (auto& y : adj_T0[x]) if (present(N[u], y) && !present(B[v], y))
						    B[w].insert(x);
					}
					for (auto& x : B[v]) {
						for (auto& y : adj_T0[x]) if (present(N[v], y) && !present(B[u], y))
								B[w].insert(x);
					}

//					db(u, v, w);
//					db(B[u], B[v], B[w]);
//					db(N[u], N[v], N[w]);

					assert(sz(B[w]) <= 2);


					si Buv;
					for (auto& x : B[u]) Buv.insert(x);
					for (auto& x : B[v]) Buv.insert(x);

					for (auto& x : Buv) TDash_bags[w].pb(x);


					assert(adj_new[u].empty() && adj_new[v].empty());
				}
			}
		}

		T.pb(mp(adj_new, root_new));
	}


//	db(T.back());
//	db(par_TDash);
	assert(n == sz(par_TDash));
	assert(T.back().sc == n - 1);



//	db(sz(adj_T0), sz(T));
//	db(sz(adj_T0), sz(TDash));

	V<si> TDashAdj(n);


	for (int i = 0; i < n; ++i) {
		if (par_TDash[i] == -1) {
			assert(TDashRoot == -1);
			TDashRoot = i;
		} else {
			add(TDashAdj, i, par_TDash[i]);
		}
	}

	assert(TDashRoot == n - 1);


	return mp(mp(TDashAdj, TDashRoot), TDash_bags);


}

void balanceTWD(VV<int>& bags, int& root, vi& par) {

	// adj_T = tree of treewidth-decomp of a CFG
	int n_T = sz(par);

	V<si> adj_T(n_T);

	function<void(int, int)> dfs = [&](int u, int parent) {
		par[u] = parent;
		for (auto& v : adj_T[u]) if (v != parent) {
				dfs(v, u);
			}
	};



	for (int i = 0; i < n_T; ++i)
		if (par[i] != -1)
			add(adj_T, i, par[i]);

	// making the tree of the treewidth-decomp binary
	{
		for (int u = 0; u < n_T; ++u) {
			if (sz(adj_T[u]) >= 4) {
				vi N;
				for (auto& elem : adj_T[u]) {
				    N.pb(elem);
				}

				int prev = u;
				for (int i = 2; i < sz(N) - 1; ++i) {
					int cur = n_T++;
					bags.pb(bags[u]);
					adj_T.pb(si());
					add(adj_T, prev, cur);
					if (i == sz(N) - 2) {
						rmv(adj_T, u, N[i]);
						rmv(adj_T, u, N[i + 1]);

						add(adj_T, cur, N[i]);
						add(adj_T, cur, N[i + 1]);
					} else {
						rmv(adj_T, u, N[i]);
						add(adj_T, cur, N[i]);
					}
					prev = cur;
				}
			}
		}
		root = -1;
		for (int i = 0; i < n_T; ++i) {
			assert(sz(adj_T[i]) <= 3);
			if (sz(adj_T[i]) <= 2)
				root = i;
		}
		assert(root != -1);

		par.assign(n_T, -1);
		dfs(root, -1);

		int minusOnes = 0;
		for (int i = 0; i < n_T; ++i) {
			if (par[i] == -1)
				++minusOnes;
		}
		assert(minusOnes == 1);
	}

	auto pr = getTWDofT(adj_T, root);


	root = pr.fs.sc;
	adj_T = pr.fs.fs;
	n_T = sz(adj_T);

	vvi bagsNew(n_T);

	for (int i = 0; i < n_T; ++i) {
		assert(sz(pr.sc[i]) <= 4);
		si bag;
		for (auto& bagId : pr.sc[i]) {
		    for (auto& u : bags[bagId]) {
		        bag.insert(u);
		    }
		}
		for (auto& elem : bag) {
			bagsNew[i].pb(elem);
		}
	}

	bags = bagsNew;

	par.assign(n_T, -1);
	dfs(root, -1);

}


int main(int arc, char *argv[]) {
#ifdef LOCAL
	auto stTime = clock();
//	freopen("../dacapo.antlr.Main.txt", "r", stdin);
//	freopen("../out.txt", "wb", stdout);
#endif
	ios::sync_with_stdio(false);
	cout.precision(10);
	cin.tie(0);
	vi badProcCntV;
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

	for (auto file: fileNames) {
		db(file);
		// take some relevant info that are independent of the analysis
		ifstream in(path + "prep_output/reachability/" + file);

		int n_G, m_G, n_H;
		vector<pair<int, int>> edgeListG;
		vector<int> vertexTypeG, procOf;
		vector<vector<vector<int>>> TWD;
		vector<int> TWD_root;
		vector<vector<int>> TWD_par;

		in >> n_G >> m_G >> n_H;


		procOf.assign(n_G, 0);

		for (int i = 0; i < n_G; ++i)
			in >> procOf[i];

		vertexTypeG.assign(n_G, 0);

		for (int i = 0; i < n_G; ++i)
			in >> vertexTypeG[i];

		edgeListG.assign(m_G, mp(-1, -1));

		for (int i = 0; i < m_G; ++i) {
			in >> edgeListG[i].fs >> edgeListG[i].sc;
			assert(edgeListG[i].fs >= 0);
			assert(edgeListG[i].fs < n_G);
			assert(edgeListG[i].sc >= 0);
			assert(edgeListG[i].sc < n_G);
		}

		;

		{
			// take TWD from pace solver
			TWD.assign(n_H, vvi());
			TWD_par.assign(n_H, vi());
			TWD_root.assign(n_H, -1);
			ifstream in(path + "treewidth_solver_output/" + file);
			string s, dummy;

			int n_BB, treewidth;
			vvi bags;
			V<si> adjTw;
			while (getline(in, s)) {
				if (s[0] == 'c')
					continue;
				stringstream ss(s);
				if (s[0] == 's') {
					ss >> dummy >> dummy;
					int nnn;
					ss >> n_BB >> treewidth >> nnn;
					--treewidth;
					assert(nnn == n_G + 1);
					adjTw.assign(n_BB, si());
				} else if (s[0] == 'b') {
					bags.pb(vi());
					ss >> dummy;
					int id;
					ss >> id;
					assert(id == sz(bags));
					int u;
					while (ss >> u) {
						--u;
						if (u != n_G) {
							assert(u >= 0);
							assert(u < n_G);
							bags.back().pb(u);
						}
					}
					// assert(!bags.back().empty());
				} else {
					int u, v;
					ss >> u >> v;
					--u, --v;
					assert(u >= 0 && v < n_BB);
					if (u == v)
						continue;
					adjTw[u].insert(v);
					adjTw[v].insert(u);
				}
			}

			{
				int cnt = 0;
				function<void(int, int)> dfs = [&](int u, int par) {
					++cnt;
					for (auto& v : adjTw[u]) if (v != par) {
							dfs(v, u);
						}
				};

				dfs(0, -1);

				assert(cnt == n_BB);
				db("adjTw is a tree!");
			}


			V<map<int, int>> bagsMapping(n_H);

			for (int i = 0; i < n_BB; ++i) {
				map<int, vi> splitting;
				for (auto& u : bags[i]) {
					assert(u >= 0 && u < n_G);
					splitting[procOf[u]].pb(u);
				}
				for (auto& pr : splitting) {
					int p = pr.fs;
					int id = sz(TWD[p]);
					bagsMapping[p][i] = id;
					TWD[p].pb(pr.sc);
				}
			}


			V<V<si>> adjTwPerProc(n_H);

			for (int p = 0; p < n_H; ++p) {
				adjTwPerProc[p].assign(sz(TWD[p]), si());
			}


			for (int i = 0; i < n_BB; ++i) {
				for (auto& j : adjTw[i]) {
					for (auto& u : bags[i]) {
						for (auto& v : bags[j]) {
							if (procOf[u] == procOf[v]) {
								int p = procOf[u];
								assert(i != j);
								assert(present(bagsMapping[p], i));
								assert(present(bagsMapping[p], j));
								assert(bagsMapping[p][i] != bagsMapping[p][j]);
								assert(bagsMapping[p][i] >= 0 && bagsMapping[p][i] < sz(TWD[p]));
								assert(bagsMapping[p][j] >= 0 && bagsMapping[p][j] < sz(TWD[p]));
								adjTwPerProc[p][bagsMapping[p][i]].insert(bagsMapping[p][j]);
								adjTwPerProc[p][bagsMapping[p][j]].insert(bagsMapping[p][i]);
							}
						}
					}
				}
			}



			db(treewidth);



			for (int p = 0; p < n_H; ++p) {
				TWD_par[p].assign(sz(TWD[p]), -1);

				function<void(int, int)> dfs = [&](int u, int par) {
					assert(u != par);
					TWD_par[p][u] = par;
					for (auto& v : adjTwPerProc[p][u]) if (v != par) {
							assert(v != u);
							assert(v != par);
							dfs(v, u);
						}
				};


				dfs(0, -1);

				int minusOnes = 0;


				for (int i = 0; i < sz(TWD[p]); ++i) {
					if (TWD_par[p][i] == -1) {
						minusOnes++;
						TWD_root[p] = i;
					}
				}
				assert(minusOnes == 1);
			}

			db("got the PACE tw-decomp");
		}

		validateTWD(TWD, TWD_root, TWD_par, procOf, edgeListG);

		vi procSz(n_H);


		for (int p = 0; p < n_H; ++p) {
			balanceTWD(TWD[p], TWD_root[p], TWD_par[p]);
		}


		validateTWD(TWD, TWD_root, TWD_par, procOf, edgeListG);


		ofstream out(path + "balanced_twd/" + file);
            int biggestBagSz = 0;
		for (int p = 0; p < n_H; ++p) {
			out << "twd " << p << '\n';
			out << sz(TWD[p]) << '\n';
			for (auto bag : TWD[p]) {
                        biggestBagSz = max(biggestBagSz, sz(bag));
				out << sz(bag);
				for (int u : bag) {
					out << ' ' << u;
					assert(procOf[u] == p);
				}
				out << '\n';
			}
			// Print twd_par
			for (int i = 0; i < sz(TWD[p]); ++i)
				out << TWD_par[p][i] << " \n"[i + 1 == sz(TWD[p])];
		}

            cout << "Balanced treewidth = " << biggestBagSz << endl;
	}


	return 0;
}
