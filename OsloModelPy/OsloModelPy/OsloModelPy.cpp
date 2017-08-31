
#include "stdafx.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <iterator>
#include <sstream>

using namespace std;

int randNum(float p) {
	float random = ((float)rand()) / (float)RAND_MAX;
	if (random <= p) {
		return 1;
	}
	else {
		return 2;
	}
}

void randsSet(int L, int *rands, float p) {
	for (int i = 0; i < L; i++) {
		rands[i] = randNum(p);
		vector<int> second(4, 100);
	}
}

extern "C" __declspec(dllexport) void iterations(int N, int L, float p, int *system, int *rands, int *avalancheSize, int *dropSize, int *heights, int *critical) {
	
	srand(time(0));
	randsSet(L, rands, p);
	bool foundCritical = false;

	for (int n = 0; n < N; n++) {
		system[0] += 1;
		bool relaxed = false;
		int s = 0;
		int d = 0;
		while (!relaxed) {
			relaxed = true;
			for (int i = 0; i < L; i++) {
				int z = system[i] - system[i + 1];
				if (z > rands[i]) {
					system[i] -= 1;
					if (i < L - 1) {
						system[i + 1] += 1;
					}
					else {
						if (!foundCritical) {
							critical[0] = n;
							foundCritical = true;
						}
						d += 1;
					}
					rands[i] = randNum(p);
					s += 1;
					relaxed = false;
				}
			}
		}
		avalancheSize[n] = s;
		dropSize[n] = d;
		heights[n] = system[0];
	}
}