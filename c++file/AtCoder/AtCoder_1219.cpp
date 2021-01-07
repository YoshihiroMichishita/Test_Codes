#include <cstdlib>
#include<iostream>
#include<string>
#include<cmath>
#define rep(i, n) for (int i = 0; i < (int)(n); i++)

using namespace std;
 
void abbreviate_fox(){
  
  int N;
  cin >> N;
  string s,t;
  cin >> s;
  rep(i,N){
    t.push_back(s[i]);
    int size = t.size();
    if (size > 2){
      if (t[size-3]=='f' && t[size-2]=='o' && t[size-1]=='x'){
        t.pop_back();
        t.pop_back();
        t.pop_back();
      }
    }
    
  }
  cout << t.size() << endl;
}

int main(int argc, char* argv[]){
  abbreviate_fox();
  return 0;
}