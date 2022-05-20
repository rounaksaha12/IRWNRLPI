#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#define ll long long

using namespace std;

ll isContains(string& str){
    
    for(ll i=0;i<(ll)str.size()-2;++i){
        if(str[i]!='>') continue;
        if(str[i+1]=='>' && str[i+2]=='>'){
            return i+2;
        }
    }
    
    return -1;
}

// read the input from 'ssearch36_rna.txt'('ssearch36_prot.txt)
int main(){
    string current;
    ll n;
    cin>>n;
    vector< vector<ll> > matrix(n, vector<ll>(n,-1));
    while(true){
        string s;
        cin>>s;
        int pos;
        if((pos=isContains(s))!=-1){
            current=string(s,pos+1);
            break;
        }
    }
    while(current!="END"){
        while(true){
            string s;
            getline(cin,s);
            if(s=="The best scores are:                                      s-w bits E(-1)"){
                break;
            }
        }

        while(true){
            string name;
            cin>>name;
            int pos;
            if((pos=isContains(name))!=-1){
                current=string(name,pos+1);
                break;
            }

            string buffer;
            char c;
            double d;
            int sw,cnt;
            string temp;

            cin>>buffer>>c>>cnt>>c>>temp>>sw>>d>>d; // for rna sequences
            // need not read temp(string) for parsing ssearch36 result of protein sequences

            ll i=stoi(current);
            ll j=stoi(name);

            cout<<i<<"\t"<<j<<"\t"<<sw<<endl;
            if(sw>matrix[i][j]){
                matrix[i][j]=sw;
            }
        }
    }

    for(ll i=0;i<matrix.size();++i){
        for(ll j=0;j<matrix[i].size();++j){
            if(matrix[i][j]==-1){
                cout<<"error"<<endl;
            }
        }
    }

    for(ll i=0;i<matrix.size();++i){
        for(ll j=i;j<matrix.size();++j){
            if(matrix[i][j]!=matrix[j][i]){
                cout<<"Asymmetric"<<endl;
                break;
            }
        }
    }

    ofstream os("rna_similarity.txt");
    // ("prot_similarity.txt")

    for(ll i=0;i<matrix.size();++i){
        for(ll j=0;j<matrix[i].size();++j){
            ll m=max(matrix[i][i],matrix[j][j]);
            os<<i<<"\t"<<j<<"\t"<<(double)matrix[i][j]/m<<endl;
        }
    }

    os.close();

    
    cout<<"end"<<endl;
    return 0;
}
