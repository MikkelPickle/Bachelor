#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace std::chrono;

random_device dev;
mt19937_64 rng(dev());
uniform_int_distribution<uint64_t> randomNumber(0,UINT64_MAX);

int64_t absoluteValue(int64_t v) {
    if (v > 0) {
        return v;
    }
    else {
        return -v;
    }
}


int nearestPowerOfTwo(int m) {
    int power = 1;
    while(power < m)
        power *= 2;
    return power;
}


class CountMinSketch {

private:

    // hash function
    // x is the input that has to be hashed
    // l is the number of counters
    // A is the random seed that has to be odd

    static uint64_t hashFunc(uint64_t x, uint64_t l, uint64_t a) {
        return (a*x) >> (64-l);
    }

public:
    int numberOfHashFunctions;
    int numberOfCounters;
    int l1norm;
    uint64_t L;
    uint64_t** grid;
    uint64_t* hashValues;

    CountMinSketch(float delta, float epsilon) {
        numberOfHashFunctions = ceil(log2(1/delta));
        l1norm = 0;
        numberOfCounters = nearestPowerOfTwo(ceil(2/epsilon));
        L = log2(numberOfCounters);
        hashValues = new uint64_t[numberOfHashFunctions];
        for (int i = 0; i < numberOfHashFunctions; i++) {
            uint64_t a = randomNumber(rng);
            a += (a % 2==0) ? 1 : 0;
            hashValues[i] = a;
        }
        grid = new uint64_t*[numberOfHashFunctions];
        for (int h = 0; h < numberOfHashFunctions; h++) {
            grid[h] = new uint64_t[numberOfCounters];
        }
    }

    void update(uint64_t item) {
        l1norm++;
        for (int i = 0; i < numberOfHashFunctions; i++) {
            uint64_t a = hashValues[i];
            uint64_t index = hashFunc(item, L, a);
            int indexToHash = static_cast<int>(index);
            grid[i][indexToHash]++;
        }
    }

    uint64_t estimate(uint64_t item) const {
        uint64_t a, index, countEstimate, newEstimate;
        a = hashValues[0];
        index = hashFunc(item, L, a);
        countEstimate = grid[0][index];
        for (int i = 1; i < numberOfHashFunctions; i++) {
            a = hashValues[i];
            index = hashFunc(item, L, a);
            newEstimate = grid[i][index];
            countEstimate = min(countEstimate, newEstimate);
        }
        return countEstimate;
    }

    void print() const {
        for (int i = 0; i < numberOfHashFunctions; i++) {
            cout << endl;
            for (int j = 0; j < numberOfCounters; j++) {
                cout << grid[i][j] << " ";
            }
        }
        cout << endl;
    }

    ~CountMinSketch() {
        for (int h = 0; h < numberOfHashFunctions; h++)
            delete[] grid[h];
        delete[] grid;
        delete[] hashValues;
    }

};



class CountSketch {

private:

    // hash function
    // x is the input that has to be hashed
    // l is the number of counters
    // A is the random seed that has to be odd
    static uint64_t hashFunc(uint64_t x, uint64_t l, uint64_t a) {
        return (a*x) >> (64-l);
    }

    static int g(uint32_t x, uint64_t a, uint64_t b) {
        uint32_t val = (a*x+b) >> 63;
        int hash = static_cast<int>(val);
        return (2*hash)-1;
    }

    static float median(int64_t* arr, int n){
        int size = n / sizeof(arr[0]);
        sort(arr, arr+size);
        if (size % 2 != 0)
            return arr[size/2];
        return (arr[(size-1)/2] + arr[size/2])/2.0;
    }

public:
    int numberOfHashFunctions;
    int numberOfCounters;
    uint64_t L;
    int64_t** grid;
    uint64_t* hashValues;
    uint64_t* aHashValues;
    uint64_t* bHashValues;

    CountSketch(float delta, float epsilon) {
        numberOfHashFunctions = ceil(log2(1/delta));
        numberOfCounters = nearestPowerOfTwo(ceil(2/epsilon));
        L = log2(numberOfCounters);
        hashValues = new uint64_t[numberOfHashFunctions];
        aHashValues = new uint64_t[numberOfHashFunctions];
        bHashValues = new uint64_t[numberOfHashFunctions];
        for (int i = 0; i < numberOfHashFunctions; i++) {
            uint64_t num = randomNumber(rng);
            num += (num % 2==0) ? 1 : 0;
            uint64_t a = randomNumber(rng);
            uint64_t b = randomNumber(rng);
            hashValues[i] = num;
            aHashValues[i] = a;
            bHashValues[i] = b;
        }
        grid = new int64_t*[numberOfHashFunctions];
        for (int h = 0; h < numberOfHashFunctions; h++) {
            grid[h] = new int64_t[numberOfCounters];
        }
    }

    void update(uint64_t item) {
        for (int i = 0; i < numberOfHashFunctions; i++) {
            uint64_t a = hashValues[i];
            uint64_t index = hashFunc(item, L, a);
            int indexToHash = static_cast<int>(index);
            auto value = static_cast<uint32_t>(item);
            uint64_t seed1 = aHashValues[i];
            uint64_t seed2 = bHashValues[i];
            int64_t count = g(value, seed1, seed2);
            grid[i][indexToHash] = grid[i][indexToHash] + count;
        }
    }

    int64_t estimate(uint64_t item) {
        int64_t a, index, estimate;
        int64_t arr [numberOfHashFunctions];
        for (int i = 0; i < numberOfHashFunctions; i++) {
            a = hashValues[i];
            index = hashFunc(item, L, a);
            estimate = grid[i][index];
            arr[i] = estimate;
        }
        return median(arr, sizeof(arr));
    }

    void print() {
        for (int i = 0; i < numberOfHashFunctions; i++) {
            cout << endl;
            for (int j = 0; j < numberOfCounters; j++) {
                cout << grid[i][j] << " ";
            }
        }
        cout << endl;
    }

    ~CountSketch() {
        for (int h = 0; h < numberOfHashFunctions; h++)
            delete[] grid[h];
        delete[] grid;
        delete[] hashValues;
    }

};




class HeavyHittersTree {

private:
    void recursiveOutput(vector<uint64_t> &heavyHitters, uint64_t v, int i) {
        uint64_t l, r;
        if (i == 0) {
            heavyHitters.push_back(v);
        }
        else {
            l = 2*v;
            r = 2*v + 1;
            if (binaryTree[i-1]->estimate(l) >= floor(epsilon * l1norm)) {
                recursiveOutput(heavyHitters,l,i-1);
            }
            if (binaryTree[i-1]->estimate(r) >= floor(epsilon * l1norm)) {
                recursiveOutput(heavyHitters,r,i-1);
            }
        }
    }

public:
    int numberOfLayersInTree;
    int l1norm;
    float epsilon;
    CountMinSketch** binaryTree;

    HeavyHittersTree(uint64_t sizeOfUniverse, float epsilon) {
        this->epsilon = epsilon;
        numberOfLayersInTree = log2(nearestPowerOfTwo(sizeOfUniverse))+1;
        binaryTree = new CountMinSketch*[numberOfLayersInTree];
        l1norm = 0;
        float epsilonPrime = epsilon/2;
        float delta = 1/pow(sizeOfUniverse, 2);
        binaryTree[0] = new CountMinSketch(delta, epsilonPrime);
        for (int i = 1; i < numberOfLayersInTree; i++) {
            binaryTree[i] = new CountMinSketch(0.25, epsilonPrime);
        }
    }

    void update(uint64_t item) {
        l1norm++;
        for (int i = 0; i < numberOfLayersInTree; i++) {
            uint64_t value = item >> i;
            binaryTree[i]->update(value);
        }
    }

    void computeHeavyHitter(vector<uint64_t> &heavyHitters) {
        cout << "Heavy Hitter threshold is " << floor(epsilon * l1norm) << endl;
        recursiveOutput(heavyHitters, 0, numberOfLayersInTree);
    }

    ~HeavyHittersTree() {
        delete[] binaryTree;
    }
};






int main() {

    /*
    vector<string> words;
    string line, tweet;

    fstream file ("cmake-build-debug/Tweets.csv", ios::in);
    if(file.is_open())
    {
        while(getline(file, line))
        {
            stringstream str(line);
            int count = 0;
            while(getline(str, tweet, ','))
                count++;
            if (count == 6) { //get column with the tweets themselves
                string word;
                stringstream str(tweet);
                while (getline(str, word, ' ')) {
                    word.erase(remove(word.begin(), word.end(), '"'), word.end());
                    word.erase(remove(word.begin(), word.end(), '!'), word.end());
                    word.erase(remove(word.begin(), word.end(), '.'), word.end());
                    if (!word.empty()) { words.push_back(word); }
                }
            }
        }
    }
    else
        cout<<"Could not open the file\n";

      //assign IDs
    unordered_map<string, int> frequencyMap;
    unordered_map<int, string> reverseMap;
    int counter = 0;
    for (auto str: words)
    {
        auto it = frequencyMap.find(str);
        if (it == frequencyMap.end()) {
            frequencyMap.insert(make_pair(str, counter));
            reverseMap.insert(make_pair(counter, str));
            counter++;
        }
    }

    //create data stream of integers that correspond to words
    vector<int> dataStream;
    for (int i = 0; i < words.size(); ++i) {
        int value = frequencyMap[words[i]];
        dataStream.push_back(value);
    }

    HeavyHittersTree h = HeavyHittersTree(words.size(), 0.01);
    for (uint64_t elem : dataStream) {
        h.update(elem);
    }
    vector<uint64_t> heavyHitters;
    h.computeHeavyHitter(heavyHitters);
    for (auto x : heavyHitters) {
        cout << reverseMap[x] << " is heavy!" << endl;
    }

     */

    vector<string> addresses;
    string line, ip;

    ifstream file ("cmake-build-debug/Darknet.CSV");
    while(getline(file, line)) {
        stringstream str(line);
        int count = 0;
        while(getline(str, ip, ',')) {
            count++;
            if (count == 4) {
                addresses.push_back(ip);
                cout << ip << endl;
            }
        }
        }
    file.close();

    //assign IDs
    unordered_map<string, int> frequencyMap;
    unordered_map<int, string> reverseMap;
    int counter = 0;
    for (auto str: addresses)
    {
        auto it = frequencyMap.find(str);
        if (it == frequencyMap.end()) {
            frequencyMap.insert(make_pair(str, counter));
            reverseMap.insert(make_pair(counter, str));
            counter++;
        }
    }

    //create data stream of integers that correspond to words
    vector<int> dataStream;
    for (int i = 0; i < addresses.size(); ++i) {
        int value = frequencyMap[addresses[i]];
        dataStream.push_back(value);
    }

    HeavyHittersTree h = HeavyHittersTree(addresses.size(), 0.01);
    for (uint64_t elem : dataStream) {
        h.update(elem);
    }
    vector<uint64_t> heavyHitters;
    h.computeHeavyHitter(heavyHitters);
    for (auto x : heavyHitters) {
        cout << reverseMap[x] << " is heavy!" << endl;
    }





    /*
    //count the precise number of occurrences
    unordered_map<int, int64_t> deterministicMap;
    CountMinSketch cms = CountMinSketch(0.10,0.01);
    for (auto i : dataStream) {
        cms.update(i);
        auto it = deterministicMap.find(i);
        if (it == deterministicMap.end()) {
            deterministicMap.insert(make_pair(i, 0));
        }
        else {
            deterministicMap.at(i) += 1;
        }
    }


    //compute the differences
    vector<int64_t> differences;
    for (auto i : deterministicMap) {
        auto estimate = cms.estimate(i.first);
        differences.push_back(absoluteValue(estimate-i.second));
    }

    sort(differences.begin(), differences.end());

    fstream file2;
    file2.open("vector_file.txt",ios_base::out);

    for(int i=0; i<differences.size(); i++)
    {
        file2<<differences[i]<<endl;
    }

    file2.close();
     */

    return 0;
}