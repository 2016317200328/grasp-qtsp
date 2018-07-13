#include <bits/stdc++.h>
#include <unistd.h>
#include <algorithm>
#include <csignal>
#include <cmath>
#include <ctime>
#include <random>

using namespace std;

#define FORR(i, a, b) for (int i = a; i < b; ++i)
#define FOR(i, a) FORR(i, 0, a)
#define MAX_SIZE 99999;
#define PROB_VALUE 1000000;

typedef pair<int,int> PairInt;
typedef vector<int> VectorInt;
typedef vector<double> VectorDouble;
typedef vector<VectorInt> VectorVI;
typedef vector<VectorDouble> VectorVD;
typedef vector<PairInt> VectorPair;

class Route {
    public:
        VectorInt route;
        double cost;
        double quota;
        Route(double c);
};

Route::Route(double c) {
    cost = c;
    quota = 0;
}

void printGRASPSolution(Route sol, VectorInt bonus) {
    cout << "=================================" << endl;
    cout << "GRASP route: ";
    
    double quota = 0.0;
    FOR(i, sol.route.size()) {
        cout << sol.route[i] << " ";
        quota += bonus[sol.route[i]];
    }
    cout << endl;

    cout << "Cost: " << sol.cost << endl;
    cout << "Quota: " << quota << endl;
    cout << "=================================" << endl << endl;
}

void printRouteIteration(VectorInt route, double cost, int ite) {
    cout << endl << "Iteration " << ite << endl;
    cout << "Route: ";
    FOR(i, route.size()) {
        cout << route[i] << " ";
    }
    cout << endl;
    cout << "Cost: " << cost << endl << endl;
}

void printRouteTwoOpt(VectorInt route) {
    cout << "Route TwoOpt: ";
    FOR(i, route.size()) {
        cout << route[i] << " ";
    }
    cout << endl << endl;
}

void printRouteSemiGreedy(VectorInt route) {
    cout << "Route SemiGreedy: ";
    FOR(i, route.size()) {
        cout << route[i] << " ";
    }
    cout << endl << endl;
}

void printProbability(VectorDouble prob) {
    cout << "Probability: ";
    FOR(i, prob.size()) {
        cout << prob[i] << "  ";
    }
    cout << endl;
}

double randDouble(double begin, double end) {
	unsigned seed = static_cast<int> (chrono::system_clock::now().time_since_epoch().count());
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(begin, end);
	return distribution(generator);
}

int compareVertex(PairInt p1, PairInt p2) {
    return p1.first - p2.first;
}

bool bonusCompare(PairInt b1, PairInt b2) {
	return b1.second >= b2.second;

}

int linearSearch(double elem, VectorDouble vector) {
    int pos = -1;

    FOR(i, vector.size()) {
        if(elem < vector[i]) {
            pos = i;
        }
    }

    return pos;
}

int findVertex(PairInt pair, VectorPair bonusPair) {
    int pos = -1;

    for(int i=0; i < bonusPair.size() && pos == -1; i++) {
        if(compareVertex(pair, bonusPair[i]) == 0) {
            pos = i;
        }
    }

    return pos;
}

bool removeVertex(VectorPair &bonusPair, PairInt pair) {
    int index = findVertex(pair, bonusPair);
    if(index >= 0 && index < bonusPair.size()) {
        bonusPair.erase(bonusPair.begin() + index);
        return true;
    }
    return false;
}

double calculateTSPCost(VectorInt tour, VectorVD costMatrix) {
    double cost = 0.0;

    FOR(i, tour.size()-1) {
        cost += costMatrix[tour[i]][tour[i+1]];
    }

    cost += costMatrix[tour.back()][tour[0]];

    return cost;
}

VectorInt twoOpt(VectorInt tour, int i, int j) {
    reverse(tour.begin() + i, tour.begin() + j+1);
    return tour;
}

VectorInt localSearchTwoOpt(VectorInt route, VectorVD costMatrix) {
    VectorInt bestRoute = route;

    FOR(i, route.size()) {
        FORR(j, i+1, route.size()) {
            VectorInt routeTwoOpt = twoOpt(bestRoute, i, j);

            //printRouteTwoOpt(routeTwoOpt);

            if (calculateTSPCost(routeTwoOpt, costMatrix) < calculateTSPCost(bestRoute, costMatrix)) {
                bestRoute = routeTwoOpt;
            }
        }
    }

    return bestRoute;
}

int roullete(VectorPair bonus) {
    double sum = 0.0;
    int probValue = PROB_VALUE;

    FOR(i, bonus.size()) {
        sum += bonus[i].second;
    }

    VectorDouble prob;
    //cout << "Sum = " << sum << endl << endl;

    prob.push_back((double) (bonus[0].second/sum));
    double maxProb = prob.back();
   // cout << "Prob[0] = " << prob.back() << " | Bonus[0] = " << bonus[0].second << endl;

    FORR(i, 1, bonus.size()) {
        //cout << "Prob[i-1] = " << prob[i-1] << " | Bonus[i].second: " << bonus[i].second;
        prob.push_back((double) ((prob[i-1]+bonus[i].second)/sum));
        //cout << " | Prob[i] = " << prob.back() << endl;
    }
    cout << endl;

    printProbability(prob);

    double minProb = prob.back();

    cout << "MaxProb = " << maxProb << " | MinProb = " << minProb << endl;

    double chosen = (double) randDouble(minProb, maxProb);

    cout << "Chosen roullete: " << chosen << endl << endl;

    //int pos = binarySearch(chosen, prob);
    int pos = linearSearch(chosen, prob);
    cout << "Pos = " << pos << endl;

    return pos;
}

VectorInt semiGreedyRoute(VectorPair bonus, int quota, int source) {
    VectorInt greedyRoute;
    greedyRoute.push_back(source);

    int sumBonus = 0;

    while(sumBonus < quota) {
        int chosen = roullete(bonus);
        greedyRoute.push_back(bonus[chosen].first);
        sumBonus += bonus[chosen].second;
        removeVertex(bonus, bonus[chosen]);
    }

    return greedyRoute;
}

void grasp(int N, VectorVD costMatrix, int quota, VectorInt bonus, int source) {
    Route bestSol(99999);
    int MAX_ITE = 1000;

    VectorPair pairBonus;
    FORR(i, 1, N) {
        pairBonus.push_back(PairInt(i, bonus[i]));
    }

    sort(pairBonus.begin(), pairBonus.end(), bonusCompare);

    cout << "BONUS ORDENADO: ";
    FOR(i, pairBonus.size()) {
        cout << pairBonus[i].second << " ";
    }
    cout << endl << endl;

    FOR(i, MAX_ITE){
        Route graspRoute(99999);
        graspRoute.route = semiGreedyRoute(pairBonus, quota, source);
        printRouteSemiGreedy(graspRoute.route);

        VectorInt routeTwoOpt;
        routeTwoOpt = localSearchTwoOpt(graspRoute.route, costMatrix);
        graspRoute.route = routeTwoOpt;

        graspRoute.cost = calculateTSPCost(graspRoute.route, costMatrix);

        //printRouteIteration(graspRoute.route, graspRoute.cost, i);

        if(graspRoute.cost < bestSol.cost) {
            bestSol.route = graspRoute.route;
            bestSol.cost = graspRoute.cost;
        }
    }

    printGRASPSolution(bestSol, bonus);
}

int main(int argc, char const *argv[]) {
    FILE * unused __attribute__((unused));
	unused = freopen (argv[1], "r", stdin);

    int N;

    cin >> N;

    VectorVD costMatrix (N, VectorDouble(N));

    FOR(i, N) {
        FOR(j, N) {
            int aux;
            cin >> aux;
            costMatrix[i][j] = aux;
        }
    }

    int quota, source = 0;

    cin >> quota;

    VectorInt bonus(N);

    FOR(i, N) {
		int vertice; int bonus_value;
		cin >> vertice >> bonus_value;
		bonus[vertice] = bonus_value;
	}

    bonus[source] = 0;

    int start = clock();
    grasp(N, costMatrix, quota, bonus, source);
    int end = clock();

    cout << "Execution time: " << (end-start)/double(CLOCKS_PER_SEC) << endl;
    return 0;
}