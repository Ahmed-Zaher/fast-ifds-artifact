#include <omp.h>
#include "template.h"
#include "IfdsInstance.h"
#include "SimpleIfdsInstance.h"
#include "algo_old.h"
#include "algo_old_demand.h"
#include "algo_new.h"


void stopClock(clock_t Time, string s = "") {

	cout << "[" << s << "]" << " time taken = " <<
	     double(clock() - Time) / CLOCKS_PER_SEC << " s" << endl;
}

class RowInExcelSheet {
public:
	string analysis, benchmark;
	ll maxD, n_G, m_G, n_GExp, m_GExp, n_H, m_H, bad_proc;
	long double preprocessing_time, all_query_time;
	long double old_tot_time, old_avg_query_time;
	long double on_demand_tot_time, on_demand_avg_query_time;

	void printRow(string path) {
		path += "result.csv";
		ofstream fOut;
		fOut.open(path, ios_base::app);
		fOut << fixed;
		fOut.precision(20);
		fOut << analysis << ","
		     << benchmark << ","
		     << maxD << ","
		     << n_G << ","
		     << m_G << ","
		     << n_GExp << ","
		     << m_GExp << ","
		     << n_H << ","
		     << m_H << ","
			 << bad_proc << ","
		     << preprocessing_time << ","
		     << all_query_time << ","
		     << old_tot_time << ","
		     << old_avg_query_time << ","
		     << on_demand_tot_time << ","
		     << on_demand_avg_query_time << "\n";

		fOut.close();
	}
};

const int TIMEOUT = 6;

int main() {
#ifdef LOCAL
	auto stTime = clock();
#endif
	ios::sync_with_stdio(false);
	cout << fixed;
	cout.precision(10);
	cin.tie(0);

	V<string> fileNames;
	string path = "../";

	{
		ofstream fOut(path + "result.csv");
		fOut << "analysis, benchmark, maxD, n_G, m_G, n_GExp, m_GExp, n_H, m_H, bad_proc, ";
		fOut << "preprocessing_time, all_query_time, old_tot_time, old_avg_query_time, ";
		fOut << "on_demand_tot_time, on_demand_avg_query_time\n";
		fOut.close();
	}


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
	fileNames.pb("avrora.txt");


	db(fileNames);


	V<string> analyses;
	// analyses.pb("reachability");
	// analyses.pb("null_ptr");
	analyses.pb("uninit_var");



	for (auto& analysis : analyses) {
		for (auto file: fileNames) {

			cout << "------------------------------------------------------------------------------------\n";
			db(file, analysis);

			RowInExcelSheet row;

			row.analysis = analysis;
			row.benchmark = file.substr(0, sz(file) - 4);

			ifstream in(path + "prep_output/" + analysis + "/" + file);

			int n_G, n_H, m_G, m_H;
			vector<pair<int, int>> edgeListG, edgeListH;
			vector<int> vertexTypeG, procOf;
			vector<vector<vector<int>>> TWD;
			vector<int> TWD_root;
			vector<vector<int>> TWD_par;
			vector<int> par_tdH;

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


			in >> m_H;

			edgeListH.assign(m_H, mp(-1, -1));

			for (int i = 0; i < m_H; ++i) {
				in >> edgeListH[i].fs >> edgeListH[i].sc;
				assert(edgeListH[i].fs >= 0);
				assert(edgeListH[i].fs < n_H);
				assert(edgeListH[i].sc >= 0);
				assert(edgeListH[i].sc < n_H);
			}

			;
			// Take the balanced treewidth decomposition
			{
				ifstream in(path + "balanced_twd/" + file);
				TWD.assign(n_H, vvi());
				TWD_par.assign(n_H, vi());
				TWD_root.assign(n_H, -1);


				for (int p = 0; p < n_H; ++p) {
					string s;
					int x;
					in >> s >> x;
					assert(s == "twd" && x == p);
					int numBags;
					in >> numBags;
					TWD[p].assign(numBags, vi());
					TWD_par[p].assign(numBags, -1);
					for (int i = 0; i < numBags; ++i) {
						int bagSz;
						in >> bagSz;
						TWD[p][i].assign(bagSz, 0);
						for (int j = 0; j < bagSz; ++j) {
							in >> TWD[p][i][j];
						}
					}
					for (int j = 0; j < numBags; ++j) {
						in >> TWD_par[p][j];
						assert(TWD_par[p][j] != j);
						assert(TWD_par[p][j] < numBags);
						if (TWD_par[p][j] == -1)
							assert(TWD_root[p] == -1), TWD_root[p] = j;
					}
				}
			}

			int treedepth;

			{

				ifstream in(path + "treedepth_solver_output/" + file);
				in >> treedepth;

				par_tdH.assign(n_H, 0);
				for (int i = 0; i < n_H; ++i) {
					if (par_tdH[i])
						assert(par_tdH[i] >= 1 && par_tdH[i] <= n_H);
					in >> par_tdH[i];
					--par_tdH[i];
				}
				vi roots;
				for (int i = 0; i < n_H; ++i) {
					if (par_tdH[i] == -1)
						roots.pb(i);
				}
				assert(sz(roots) == 1);
			}

			db(treedepth);

			vi D(n_H);

			for (int i = 0; i < n_H; ++i)
				in >> D[i];


			int m_GExp = 0;
			vector<vector<vector<int>>> flow(m_G);

			for (int e = 0; e < m_G; ++e) {
				int u = edgeListG[e].fs, v = edgeListG[e].sc;
				assert(u < n_G && u < sz(procOf));
				assert(procOf[u] < n_H);
				flow[e].assign(D[procOf[u]] + 1, vi());

				for (int j = 0; j <= D[procOf[u]]; ++j) {
					int Sz;
					in >> Sz;

					for (int k = 0; k < Sz; ++k) {
						int jj;
						in >> jj;
						assert(jj >= 0 && jj <= D[procOf[v]]);
						flow[e][j].pb(jj);
						++m_GExp;
					}
				}
			}

			string done;
			in >> done;

			assert(done == "done");



			assert(m_G == sz(edgeListG));
			assert(m_H == sz(edgeListH));
			clock_t t = clock();

			IfdsInstance instance(n_G, n_H, m_G, m_H, D, edgeListG, edgeListH, flow, vertexTypeG, procOf, TWD,
			                      TWD_root, TWD_par, par_tdH);


			row.maxD = 0;

			for (int i = 0; i < instance.n_H; ++i) {
				row.maxD = max(row.maxD, 1LL * instance.D[i]);
			}

			row.n_G = instance.n_G;
			row.m_G = instance.m_G;
			row.n_GExp = instance.n_GExp;
			row.m_GExp = instance.m_GExp;
			row.n_H = instance.n_H;
			row.m_H = instance.m_H;
			row.bad_proc = instance.badProcCnt;


			stopClock(t, "Creating instance time");



			db(instance.badProcCnt);
			db(instance.badProcPercentage);


			t = clock();
			algo_new alg(instance);
			clock_t totClockCyclesPreprocessing = clock() - t;
			stopClock(t, "Total preprocessing time");


			row.preprocessing_time = ((long double) totClockCyclesPreprocessing) / CLOCKS_PER_SEC;

			{
				vi returnSite(instance.n_G, -1);
				for (int i = 0; i < instance.n_G; ++i)
					if (instance.vertexTypeG[i] == instance.CALL_VERTEX) {
						assert(sz(instance.G[i]) == 2);
						for (auto& v : instance.G[i]) {
							if (instance.vertexTypeG[v] == instance.RETURN_SITE_VERTEX)
								assert(returnSite[i] == -1), returnSite[i] = v;
						}
						assert(returnSite[i] != -1);
						assert(instance.procOf[i] == instance.procOf[returnSite[i]]);
					}
			}

			cout << "------------------------------------------------------------------------------------\n";


			map<int, vi> querySet = instance.generateQuerySet();
			clock_t totClockCyclesOld = 0;
			clock_t totClockCyclesOldDemand = 0;
			clock_t totClockCyclesAlg = 0;
			long double totBytesUsedOld = 0;
			long double totBytesUsedDemand = 0;
			int queryCntOld = 0;
			int queryCntDemand = 0;


			// passing queries to our algo
			for (int U = 0; U < instance.n_GExp; ++U) {
				for (auto& V : querySet[U]) {
					t = clock();
					bool ans = alg.query(U, V);
					totClockCyclesAlg += clock() - t;
				}
			}

			row.all_query_time = ((long double) totClockCyclesAlg) / CLOCKS_PER_SEC;

			SimpleIfdsInstance simpleInstance = SimpleIfdsInstance(instance);

			// passing queries to old algo

			for (int U = 0; U < instance.n_GExp; ++U)
				if (!querySet[U].empty()) {
					for (auto& V : querySet[U]) {
						++queryCntOld;
						algo_old alg_oldie(&simpleInstance, U);
						totClockCyclesOld += alg_oldie.runTime;
						t = clock();
						bool ansOld = alg_oldie.query(V);
						totClockCyclesOld += clock() - t;
						db(queryCntOld, double(alg_oldie.runTime) / CLOCKS_PER_SEC);
						db(queryCntOld, double(totClockCyclesOld) / CLOCKS_PER_SEC);
						if (((long double) totClockCyclesOld) / CLOCKS_PER_SEC > TIMEOUT)
							break;
					}
					if (((long double) totClockCyclesOld) / CLOCKS_PER_SEC > TIMEOUT)
						break;
				}


			row.old_tot_time = ((long double) totClockCyclesOld) / CLOCKS_PER_SEC;
			row.old_avg_query_time = ((long double) totClockCyclesOld) / CLOCKS_PER_SEC / queryCntOld;




			// passing queries to demand version of old algo

			for (int U = 0; U < instance.n_GExp; ++U)
				if (!querySet[U].empty()) {
					for (auto& V : querySet[U]) {
						++queryCntDemand;
						algo_old_demand alg_oldie_demand(&simpleInstance, U);
						totClockCyclesOldDemand += alg_oldie_demand.runTime;
						t = clock();
						bool ansOld = alg_oldie_demand.query(V);
						totClockCyclesOldDemand += clock() - t;

						db(queryCntDemand, double(alg_oldie_demand.runTime) / CLOCKS_PER_SEC);
						db(queryCntDemand, double(totClockCyclesOldDemand) / CLOCKS_PER_SEC);
						if (((long double) totClockCyclesOldDemand) / CLOCKS_PER_SEC > TIMEOUT)
							break;
					}
					if (((long double) totClockCyclesOldDemand) / CLOCKS_PER_SEC > TIMEOUT)
						break;
				}



			row.on_demand_tot_time = ((long double) totClockCyclesOldDemand) / CLOCKS_PER_SEC;
			row.on_demand_avg_query_time = ((long double) totClockCyclesOldDemand) / CLOCKS_PER_SEC / queryCntDemand;






			// IMP IMP IMP
			row.printRow(path);

			cout << "------------------------------------------------------------------------------------\n";
		}

	}


#ifdef LOCAL
	cout << "\n\n\nExecution time: " <<
	     (clock() - stTime) * 1e3 / CLOCKS_PER_SEC << " ms" << endl;
#endif
	return 0;
}
