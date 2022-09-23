#ifndef CPP_CODE_SIMPLEIFDSINSTANCE_H
#define CPP_CODE_SIMPLEIFDSINSTANCE_H

#include "template.h"
#include "IfdsInstance.h"

class SimpleIfdsInstance {
public:
	int n_G, n_H, n_GExp;
	vector<vector<int>> G, GExp;
	vector<int> vertexTypeG, procOf;
	vector<int> D;
	vector<vector<int>> nodeGExp;
	vector<pair<int, int>> revMapGExp;
	vi s;

	SimpleIfdsInstance(IfdsInstance instance) {
		n_G = instance.n_G;
		n_H = instance.n_H;
		n_GExp = instance.n_GExp;
		D = instance.D;
		G = instance.G;
		GExp = instance.GExp;
		vertexTypeG = instance.vertexTypeG;
		procOf = instance.procOf;
		nodeGExp = instance.nodeGExp;
		revMapGExp = instance.revMapGExp;
		s = instance.s;
	}
};


#endif //CPP_CODE_SIMPLEIFDSINSTANCE_H
