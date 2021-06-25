#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <queue>
//#include <bits/stdc++.h>
#include <chrono>
#include <set>
#include <cassert>
#include <iomanip>
//#include <map>
#include <cmath>
#include <bitset>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <string.h>

using namespace std;

#define REP(i,x,v)for(int i=x;i<=v;i++)
#define REPD(i,x,v)for(int i=x;i>=v;i--)
#define FOR(i,v)for(int i=0;i<v;i++)

#define SZ(x) ((int)(x).size())
#define pb push_back

#define sim template < class c
#define ris return * this
#define dor > debug & operator <<
#define eni(x) sim > typename   enable_if<sizeof dud<c>(0) x 1, debug&>::type operator<<(c i) {
sim > struct rge { c b, e; };
sim > rge<c> range(c i, c j) { return rge<c>{i, j}; }
sim > auto dud(c* x) -> decltype(cerr << *x, 0);
sim > char dud(...);
struct debug {
~debug() { cerr << endl; }
eni(!=) cerr << boolalpha << i; ris; }
eni(==) ris << range(begin(i), end(i)); }
sim, class b dor(pair < b, c > d) {
  ris << "(" << d.first << ", " << d.second << ")";
}
sim dor(rge<c> d) {
  *this << "[";
  for (auto it = d.b; it != d.e; ++it)
    *this << ", " + 2 * (it == d.b) << *it;
  ris << "]";
}
};
#define var(...) " [" << #__VA_ARGS__ ": " << (__VA_ARGS__) << "] "

using vpvii=vector<pair<vector<int>, int>>;
using ll=long long;
using pii=pair<int,int>;
using pll=pair<ll,ll>;
using vi=vector<int>;
using vd=vector<double>;
using vll=vector<ll>;


mt19937_64 rng(134569);
int rd(int l, int r) {return uniform_int_distribution<int>(l, r)(rng);}

#define FORESTSIZE 37
#define SUNCYCLE 6
#define INF (1<<30)
#define NGBHNO 6
#define MAXDAY 23
#define PLEN 12
#define PENLEN 24
#define TIMETOSTOP 90.0
#define MAXACTIONS 38

using namespace std;

using uns = unsigned long long;
const uns X=5111117;
using ull = unsigned long long;


ull allAvailableTrees;
double pen[PENLEN];
int forestMap[FORESTSIZE][NGBHNO];
int soilRichness[FORESTSIZE];
int distFromTree[FORESTSIZE][FORESTSIZE];
string actionName[4]= {"SEED", "GROW", "COMPLETE", "WAIT"};
int basePrices[5] = {0,1,3,7,4};
int richnessBonus[4] = {0,0,2,4};
chrono::steady_clock::time_point beginTurn;
ull treeRange[FORESTSIZE][4];
int sunPointsComputed[2] = {-1, -1};
ull getsShadowedBy[FORESTSIZE][6][4][4];

//TODO get rid of expensive asserts in state

void computePenalty(int day) {
    // double lambda = 1.0 - day/400.0;
    double lambda;
    // if(day<18) lambda = 0.98;
    // else lambda = 0.95;
    double l0 = 0.99;
    double l1 = 0.94;
    double t = (1.0*day*day)/(1.0*MAXDAY*MAXDAY);
    lambda = l0 * (1.0 - t) + t * l1;
    pen[0] = 1;
    for(int i=1; i< PENLEN; i++){
        pen[i] = pen[i-1] * lambda;
    }
}






vector<int> computeDistfromTree(int treeIdx){
    vector<int> dist(FORESTSIZE, -1);
    queue<int> que;
    que.push(treeIdx);
    dist[treeIdx] = 0;
    int cnt = 0;
    while(!que.empty()){
        cnt += 1;
        int currIdx = que.front();
        que.pop();
        for(int n=0; n<NGBHNO; n++){
            int ngbIdx = forestMap[currIdx][n];
            if(ngbIdx!=-1 && dist[ngbIdx]==-1){
                dist[ngbIdx] = dist[currIdx] + 1;
                que.push(ngbIdx);
            }
        }
    }
    return dist;
}

void precomputeDistancesFromTrees(){
    for(int i=0; i<FORESTSIZE; i++){
        vector<int> dist = computeDistfromTree(i);
        for(int j=0; j<FORESTSIZE; j++) distFromTree[i][j] = dist[j];
    }
}



#define SETBIT(x, pos) x|=(1LL<<(pos))
#define UNSETBIT(x, pos) x&=(~(1LL<<(pos)))



struct Hex{
    int* a;

    Hex(int *_a){
        a = _a;
    }
};

ostream& operator<<(ostream& os, const Hex& h) {
    os << endl;
    int printOrder[FORESTSIZE] = {25,24,23,22,26,11,10,9,21,27,12,3,2,8,20,28,13,4,0,1,7,19,29,14,5,6,18,36,30,15,16,17,35,31,32,33,34};
    int endLine[FORESTSIZE] = {0};
    int newLine[FORESTSIZE] = {0};
    int endLineIdx[7] = {3, 8, 14, 21, 27, 32, 36};
    int newLineIdx[7] = {0, 4, 9, 15, 22, 28, 33};
    for(auto ei: endLineIdx) endLine[ei] = 1;
    for(auto nl: newLineIdx) newLine[nl] = 1;
    int l=0;
    for(int i=0; i<FORESTSIZE; i++ ){
        int extraSpace =0;
        if(newLine[i]){ extraSpace = 7-(endLineIdx[l] - newLineIdx[l] + 1); l+=1;}
        string s;
        (endLine[i])? s= "\n" : s=" ";
        for(int sp=0; sp < extraSpace; sp++) {os << " "; }
        (h.a[printOrder[i]]!=-1)? os << h.a[printOrder[i]] : os << "*";
        os << s;
    }
    return os;
}



struct Action{
    int type;  // 0-seed, 1-grow, 2-complete, 3-wait
    int source;
    int target;

    Action(){
        type = -1;
        source = -1;
        target = -1;
    }

    Action(int _type, int _source, int _target){
        type = _type;
        source = _source;
        target = _target;
    }
};

Action waitAction() {
    return Action(3,-1,-1);
}


struct Tree{
    int size;
    int owner;//-1, 0, 1
    bool isDormant;

    Tree(int _size, int _owner, bool _isDormant){
        size = _size;
        owner = _owner;
        isDormant = _isDormant;
    }

    Tree(){
        size = -1;
        owner = -1;
        isDormant = false;
    }
};

ostream& operator<<(ostream& os, const Tree &t) {
        os<<"Tree: sz "<<t.size<<" owner "<<t.owner<<" dorm "<<t.isDormant;
        return os;
    }

ostream& operator<<(ostream& os, const Action &a) {
        if(a.type == 0){
            os << actionName[0] << " " << a.target << " " << a.source ;
        }
        if(a.type == 1 || a.type == 2){
            os << actionName[a.type] << " "  << a.source ;
        }
        if(a.type == 3){
            os << "WAIT" ;
        }
        return os;
    }


int treeIdxs[2][FORESTSIZE];
int treeCnts[2];

struct GameState{
    ull treeOwnedBy[2];
    ull treeSize[4];
    ull treeDormant;
    int day;
    int nutritient;
    int sunPoints[2];
    int score[2];
    unsigned isWaiting;

    GameState(){}

    GameState(int _day, int _nutritient, int _sunPoints, int _score, int _oppSunPoints, int _oppScore, bool _oppIsWaiting){
        FOR(i, 4)treeSize[i]=0;
        FOR(i, 2) treeOwnedBy[i]=0;
        treeDormant = 0;
        day = _day;
        nutritient = _nutritient;
        sunPoints[0] = _sunPoints;
        sunPoints[1] = _oppSunPoints;
        score[0] = _score;
        score[1] = _oppScore;
        isWaiting=0;
        isWaiting=2*((int)_oppIsWaiting);
    }

    inline void zeroTree(int cell) {
        FOR(p, 2) {UNSETBIT(treeOwnedBy[p], cell);}
        FOR(s, 4) {UNSETBIT(treeSize[s], cell);}
        UNSETBIT(treeDormant, cell);
    }

    inline bool isPlayerWaiting(int player) {
        return isWaiting&(1<<player);
    }

    inline bool isTreeDormant(int cell) {
        return treeDormant&(1LL<<cell);
    }

    inline int countTreesOfSize(int player, int size) {
        return __builtin_popcountll(treeOwnedBy[player]&treeSize[size]);
    }

    inline int costOfSeed(int player){
        return countTreesOfSize(player, 0);
    }

    int costOfGrowTo(int player, int size){
        assert(size>0 && size <4);
        return basePrices[size] + countTreesOfSize(player, size);
    }

    void addTree(int cell, Tree t){
        zeroTree(cell);
        SETBIT(treeOwnedBy[t.owner], cell);
        SETBIT(treeSize[t.size], cell);
        if (t.isDormant) SETBIT(treeDormant, cell);
    }

    inline int getTreeOwner(int cell) {
        if (treeOwnedBy[0]&(1LL<<cell)) return 0;
        if (treeOwnedBy[1]&(1LL<<cell)) return 1;
        return -1;
    }

    inline int getTreeSize(int cell) {
        REPD(s, 3, 0) {
            if (treeSize[s]&(1LL<<cell)) return s;
        }
        assert(false);
    }

    void growTree(int player, int cell){
        //assert(player == getTreeOwner(cell));
        int currSize = getTreeSize(cell);
        sunPoints[player] -= costOfGrowTo(player, currSize + 1);
        UNSETBIT(treeSize[currSize], cell);
        SETBIT(treeSize[currSize+1], cell);
        SETBIT(treeDormant, cell);
    }

    void seedTree(int player, int source, int target) {
        //assert(player == getTreeOwner(source));
        sunPoints[player] -= costOfSeed(player);
        SETBIT(treeDormant, source);
        SETBIT(treeOwnedBy[player], target);
        SETBIT(treeSize[0], target);
    }

    void completeTree(int player, int cell) {
        //assert(player == getTreeOwner(cell));
        //assert(soilRichness[cell] > 0);
        int gain = nutritient + richnessBonus[soilRichness[cell]];
        score[player] += gain;
        nutritient = max(0, nutritient-1);
        zeroTree(cell);
        sunPoints[player] -= 4;
    }

    void wait(int player){
        isWaiting |= (1<<player);
    }

    void newDay() {
        day += 1;  // zeby dobrze liczyl sie shadow, dzien trzeba zwiekszyc na poczatku a nie na koncu!!!!!
        isWaiting = 0;
        if(day<MAXDAY +1) {
            treeDormant = 0;
            int sunDir = (day+3)%6;
            FOR(p, 2) {
                REP(s, 1, 3) {
                    ull treeMask = treeSize[s]&treeOwnedBy[p];
                    while(treeMask) {
                        int idx = ffsll(treeMask)-1;
                        treeMask ^= (1LL<<idx);
                        bool shadowed = false;
                        REP(s2, s, 3) {
                            if (getsShadowedBy[idx][sunDir][s][s2]&treeSize[s2]) {shadowed=true;break;}
                        }
                        if(!shadowed) sunPoints[p] += s;
                    }
                }
            }

        }
    }

    bool isTerminal() {
        return day == 24;
    }

    double reward() {
        int totTrees[2] = {0, 0};
        FOR(p, 2) FOR(i, 4) totTrees[p]+=countTreesOfSize(p, i);
        return (score[0] - score[1]) + (sunPoints[0]/3 - sunPoints[1]/3) + 0.01*(totTrees[0] - totTrees[1]);
    }

    void applyAction(int player, const Action& a) {
        if (a.type == 0) seedTree(player, a.source, a.target);
        if (a.type == 1) growTree(player, a.target);
        if (a.type == 2) completeTree(player, a.target);
        if (a.type == 3) wait(player);
    }

    void applyActions(const Action& a0, const Action& a1) {
        // both seed
        if (a0.type == 0 && a1.type == 0 && a0.target == a1.target) {
            SETBIT(treeDormant, a0.source);
            SETBIT(treeDormant, a1.source);
            return;
        }

        // both complete
        if (a0.type == 2 && a1.type == 2) {
            int oldnutritient = nutritient;
            applyAction(0, a0);
            nutritient = oldnutritient;
            applyAction(1, a1);
            nutritient = max(0, oldnutritient - 2);
            return;
        }

        // both wait
        if (a0.type == 3 && a1.type == 3) {
            newDay();
            return;
        }
        applyAction(0, a0);
        applyAction(1, a1);
    }

    double rollout() {
        int sunPForWait[2][6];
        for(int i=0; i<6; i++){
            int d = day;
            if(d>=MAXDAY) break;
            int prevSP[2] = {sunPoints[0], sunPoints[1]};
            newDay();
            FOR(p,2) sunPForWait[p][d%6] = sunPoints[p] - prevSP[p];
        }

        while(day <MAXDAY) {
            FOR(p,2) sunPoints[p] += sunPForWait[p][day%6];
            day+=1;
        }

        if (day == MAXDAY) {

            FOR(p,2) treeCnts[p]=0;

            for(int i=0; i<FORESTSIZE; i++) {
                int p = getTreeOwner(i);
                if (p==-1) continue;
                if ((treeSize[3]&(1LL<<i)) && !isTreeDormant(i)) {
                    treeIdxs[p][treeCnts[p]++] = i;
                }
            }

            for(int i=0;;i++) {
                if (i>=treeCnts[0] && i>=treeCnts[1]) break;
                int cntComplete=0;
                FOR(p, 2) {
                    if (i<treeCnts[p] && sunPoints[p]>=4) {
                        int idx = treeIdxs[p][i];
                        int pointGain = nutritient + richnessBonus[soilRichness[idx]];
                        if (3*pointGain > 4) {
                            sunPoints[p] -= 4;
                            score[p] += pointGain;
                            cntComplete++;
                            zeroTree(idx);
                        }
                    }
                }
                nutritient = max(0, nutritient-cntComplete);
                if (cntComplete == 0) break;
            }
        }

        return reward();
    }

    uns hash() const {
        uns z=0;
        FOR(i, 2) z = z*X+treeOwnedBy[i];
        FOR(i, 4) z = z*X+treeSize[i];
        z = z*X + treeDormant;
        z = z*X+unsigned(day);
        z = z*X+unsigned(nutritient);
        z = z*X+unsigned(sunPoints[0]);
        z = z*X+unsigned(sunPoints[1]);
        z = z*X+unsigned(score[0]);
        z = z*X+unsigned(score[1]);
        z = z*X+isWaiting;
        return z;
    }




    void genPossibleActions(int player, Action (&aa)[MAXACTIONS], int &cnt) {
        cnt=0;
        aa[cnt++] = waitAction();

        if (isPlayerWaiting(player)) {
            return;
        }

        if (sunPoints[player]>=4) {
            ull size3 = treeSize[3]&treeOwnedBy[player]&(~treeDormant);
            while(size3) {
                int idx = ffsll(size3)-1;
                size3 ^= (1LL<<idx);
                aa[cnt++] = Action(2, idx, idx);
            }
        }

        REP(s, 0, 2) {
            if (sunPoints[player]>=costOfGrowTo(player, s+1)) {
                ull treeMask = treeSize[s]&treeOwnedBy[player]&(~treeDormant);
                while(treeMask) {
                    int idx = ffsll(treeMask)-1;
                    treeMask ^= (1LL<<idx);
                    aa[cnt++] = Action(1, idx, idx);
                }
            }
        }


        if (sunPoints[player] >= costOfSeed(player)) {
            ull seedable = allAvailableTrees&(~treeOwnedBy[0])&(~treeOwnedBy[1]);
            REPD(s, 3, 1) {
                ull treeMask = treeSize[s]&treeOwnedBy[player]&(~treeDormant);
                while(treeMask) {
                    int i = ffsll(treeMask)-1;
                    treeMask ^= (1LL<<i);
                    ull newCells = seedable & treeRange[i][s];
                    seedable ^= newCells;
                    while (newCells > 0) {
                        int idx = ffsll(newCells)-1;
                        newCells ^= (1LL<<idx);
                        aa[cnt++] = Action(0, i, idx);
                    }
                }
            }
        }
    }
};



ostream& operator<<(ostream& os, const GameState &s) {
    os << "Day: " << s.day << ", ";
    os<< "Nutritient: " << s.nutritient << ", ";
    os << "SunP: " << s.sunPoints[0]<< ", " <<s.sunPoints[1]<<" ";
    os << "Score: " << s.score[0] << ", "<< s.score[1]<<" ";
    return os;
}


struct BeamNode {
    GameState state;
    Action firstAction;
    double currValue;
    double fullValue;

    BeamNode(GameState _state, Action _firstAction, double _currValue, double _fullValue){
        state = _state;
        firstAction = _firstAction;
        currValue = _currValue;
        fullValue = _fullValue;
    }


    bool operator < (const BeamNode& beamNode2) const {
    return fullValue < beamNode2.fullValue;
    }
};



double nodeValue(const BeamNode& node) {
    double fullValue = node.currValue;
    GameState s = node.state;
    int sunPForWait[2][6];
    for(int i=0; i<6; i++){
        int d = s.day;
        if(d>=MAXDAY) break;
        int prevSP[2] = {s.sunPoints[0], s.sunPoints[1]};
        s.newDay();
        FOR(p, 2) sunPForWait[p][d%6] = s.sunPoints[p] - prevSP[p];
    }

    while(s.day <MAXDAY) {
        FOR(p,2) s.sunPoints[p] += sunPForWait[p][s.day%6];
        s.day+=1;
    }

    if (s.day == MAXDAY) {
        ull treeMask = s.treeSize[3]&s.treeOwnedBy[0]&(~s.treeDormant);
        while(treeMask) {
            if (s.sunPoints[0]<4) break;
            int idx = ffsll(treeMask)-1;
            treeMask ^= (1LL<<idx);
            int pointGain = s.nutritient + richnessBonus[soilRichness[idx]];
            if (3*pointGain > 4) {
                s.sunPoints[0] -= 4;
                s.score[0] += pointGain;
                fullValue += pointGain*pen[MAXDAY];
                s.nutritient = max(0, s.nutritient-1);
            }
        }
    }

    fullValue+=(s.sunPoints[0]-s.sunPoints[1])*pen[MAXDAY]*(0.25);
    return fullValue;

}





vector<BeamNode> adjacentNodes;
unordered_map<uns, double> nodeMem;
int cntNodes;
int cntUniqueNodes;


Action expActions[MAXACTIONS];


void expandNode(BeamNode node, priority_queue<BeamNode> &pq) {
    if(node.state.day==MAXDAY+1){
        pq.push(node);
        return;
    }
    assert(node.state.day <= MAXDAY);

    int nActions=0;
    node.state.genPossibleActions(0, expActions, nActions);
    if(node.firstAction.type == -1){
        vector<Action> va(expActions, expActions + nActions);
        debug() << var(va);
    }

    FOR(i, nActions) {
        cntNodes++;
        BeamNode newNode = node;
        const Action& a = expActions[i];
        newNode.state.applyActions(a, waitAction());
        int gain = newNode.state.score[0] - node.state.score[0];
        newNode.currValue += gain*pen[node.state.day];
        newNode.fullValue = nodeValue(newNode);
        unsigned newHash = newNode.state.hash();
        if(nodeMem.count(newHash) > 0 && nodeMem[newHash]>= newNode.fullValue ) {
            return;
        }
        cntUniqueNodes++;
        nodeMem[newHash] = newNode.fullValue;
        if (node.firstAction.type == -1) {
            newNode.firstAction = a;
        }
        pq.push(newNode);
    }
}




Action findAction(GameState s){
    computePenalty(s.day);
    nodeMem.clear();
    nodeMem.reserve(10000);
    cntNodes = 0;
    cntUniqueNodes = 0;

    Action chosenAction;
    priority_queue<BeamNode> pq[PLEN+1];
    cerr << "jestem w findAction" << endl;
    pq[0].push(BeamNode(s,Action(-1,-1,-1), 0.0, 0.0)); //
    int iter = 0;
    for(;;iter++) {
        chrono::steady_clock::time_point checkPoint = chrono::steady_clock::now();
        auto timeFromTurnStart = chrono::duration_cast<chrono::milliseconds>(checkPoint - beginTurn).count();
        if (timeFromTurnStart>TIMETOSTOP) break;
        for(int l=0; l<PLEN; l++){
            if(pq[l].empty()) continue; // chyba?
            BeamNode node = pq[l].top();
            pq[l].pop();
            expandNode(node, pq[l+1]);
        }
    }
    debug()<<var(iter);
    debug() << var(cntNodes);
    debug() << var(cntUniqueNodes);

    assert(!pq[PLEN].empty());
    chosenAction= pq[PLEN].top().firstAction;
    double val =  pq[PLEN].top().fullValue;
    debug()<<var(val);

    return chosenAction;
}

void precomputeTreeRange() {
    for(int tree=0;tree<FORESTSIZE;tree++) {
        for(int d=0;d<=3;d++) {
            ull inRange = 0;
            for(int u=0;u<FORESTSIZE;u++) {
                if (distFromTree[tree][u]<=d) {
                    inRange += (1LL<<u);
                }
            }
            treeRange[tree][d]=inRange;
        }
    }
}

void precomputeShadows() {
    for(int i=0; i<FORESTSIZE; i++) {
        FOR(sunDir, 6) {
            REP(s, 1, 3) {
                int currIdx = i;
                for(int d=1; d<=3; d++) {
                    currIdx = forestMap[currIdx][sunDir];
                    if(currIdx == -1) break;
                    REP(s2, max(d, s), 3) {
                        getsShadowedBy[i][sunDir][s][s2]|=(1LL<<currIdx);
                    }
                }
            }
        }
    }
}




void precompute(){
    precomputeDistancesFromTrees();
    precomputeTreeRange();
    precomputeShadows();
    // computePenalty();


}

void readInitialData(){
    allAvailableTrees = 0;
    int numberOfCells; // 37
    cin >> numberOfCells; cin.ignore();
    for (int i = 0; i < numberOfCells; i++) {
        int index; // 0 is the center cell, the next cells spiral outwards
        int richness; // 0 if the cell is unusable, 1-3 for usable cells
        int neigh0; // the index of the neighbouring cell for each direction
        int neigh1;
        int neigh2;
        int neigh3;
        int neigh4;
        int neigh5;
        cin >> index >> richness >> neigh0 >> neigh1 >> neigh2 >> neigh3 >> neigh4 >> neigh5; cin.ignore();
        if (richness>0) {
            allAvailableTrees |= (1LL<<i);
        }

        soilRichness[i] = richness;
        forestMap[i][0] = neigh0;
        forestMap[i][1] = neigh1;
        forestMap[i][2] = neigh2;
        forestMap[i][3] = neigh3;
        forestMap[i][4] = neigh4;
        forestMap[i][5] = neigh5;
    }
}

void play(){


    // game loop
    while (1) {


        int day; // the game lasts 24 days: 0-23
        cin >> day; cin.ignore();
        int nutrients; // the base score you gain from the next COMPLETE action
        cin >> nutrients; cin.ignore();
        int sun; // your sun points
        int score; // your current score
        cin >> sun >> score; cin.ignore();
        int oppSun; // opponent's sun points
        int oppScore; // opponent's score
        bool oppIsWaiting; // whether your opponent is asleep until the next day
        cin >> oppSun >> oppScore >> oppIsWaiting; cin.ignore();

        GameState state = GameState(day, nutrients, sun, score, oppSun, oppScore, oppIsWaiting);

        int numberOfTrees; // the current amount of trees
        cin >> numberOfTrees; cin.ignore();
        for (int i = 0; i < numberOfTrees; i++) {
            int cellIndex; // location of this tree
            int size; // size of this tree: 0-3
            bool isMine; // 1 if this is your tree
            bool isDormant; // 1 if this tree is dormant
            cin >> cellIndex >> size >> isMine >> isDormant; cin.ignore();

            int owner = isMine? 0: 1;
            state.addTree(cellIndex, Tree(size, owner, isDormant));
        }


        int numberOfPossibleActions; // all legal actions
        cin >> numberOfPossibleActions; cin.ignore();
        for (int i = 0; i < numberOfPossibleActions; i++) {
            string possibleAction;
            getline(cin, possibleAction); // try printing something from here to start with
        }



        if(sunPointsComputed[0] != -1) {
            //cerr << "SERIO TO JEST DOBRZE????????" <<endl;
            debug() << var(sunPointsComputed[0]) << var(state.sunPoints[0]);
            debug() << var(sunPointsComputed[1]) << var(state.sunPoints[1]);
            assert(sunPointsComputed[0] == state.sunPoints[0]);
            assert(sunPointsComputed[1] == state.sunPoints[1]);
        }

        GameState copyState = state;
        cerr<< "poczatek tury"<< endl;
        beginTurn = chrono::steady_clock::now();


        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        debug() << var(state);
        Action myAction = findAction(state);

        if(myAction.type == 3 && oppIsWaiting && state.day<MAXDAY){
            copyState.newDay();
            FOR(p, 2) sunPointsComputed[p] = copyState.sunPoints[p];
        }
        else {
            FOR(p, 2) sunPointsComputed[p] = -1;
        }

        if(myAction.type == 0){
            cout << actionName[0] << " " <<myAction.source << " " <<  myAction.target << endl;
        }
        if(myAction.type == 1 || myAction.type == 2){
            cout << actionName[myAction.type] << " "  << myAction.source << endl;
        }
        if(myAction.type == 3){
        cout << "WAIT" << endl;
        }
        assert(myAction.type != -1);
        assert(myAction.type != 4);


        // GROW cellIdx | SEED sourceIdx targetIdx | COMPLETE cellIdx | WAIT <message>
    }
}


int main(){
    readInitialData();
    precompute();
    play();
}