/*
 * DP algorithm to count labeled chordal graphs.
 */

#include<iostream>
#include<vector>
#include<string>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;
using namespace std;

typedef int1024_t lli;

bool verbatim = false;
int rec_depth = 0;
lli notSaved = -1;
int w;

lli choose(int n, int k);

lli chordal(int k);
lli comp_chordal(int k);
vector<lli> table_chordal;

lli chordal_conn(int k);
lli comp_chordal_conn(int k);
vector<lli> table_chordal_conn;

lli g(int t, int x, int z, int k);
lli comp_g(int t, int x, int z, int k);
vector<vector<vector<vector<lli> > > > table_g;

lli gTilde(int t, int x, int z, int k);
lli comp_gTilde(int t, int x, int z, int k);
vector<vector<vector<vector <lli> > > > table_gTilde;

lli gHat(int t, int x, int z, int k);
lli comp_gHat(int t, int x, int z, int k);
vector<vector<vector<vector <lli> > > > table_gHat;

lli g1Tilde(int t, int x, int k);
lli comp_g1Tilde(int t, int x, int k);
vector<vector<vector<lli> > > table_g1Tilde;

lli g2Tilde(int t, int x, int k);
lli comp_g2Tilde(int t, int x, int k);
vector<vector<vector<lli> > > table_g2Tilde;

lli f(int t, int x, int l, int k);
lli comp_f(int t, int x, int l, int k);
vector<vector<vector<vector <lli> > > > table_f;

lli fTilde(int t, int x, int l, int k);
lli comp_fTilde(int t, int x, int l, int k);
vector<vector<vector<vector <lli> > > > table_fTilde;

lli fHat(int t, int x, int l, int k);
lli comp_fHat(int t, int x, int l, int k);
vector<vector<vector<vector <lli> > > > table_fHat;

lli fHat(int t, int x, int z, int l, int k);
lli comp_fHat(int t, int x, int z, int l, int k);
vector<vector<vector<vector <vector<lli> > > > > table_fHat5;

vector<vector<lli> > table_choose;


// DEBUG PRINTOUTS
void debugIn(string fname, lli x1, lli x2, lli x3) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << ")" << endl;}

void debugIn(string fname, lli x1, lli x2, lli x3, lli x4) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << "," << x4 << ")" << endl;}

void debugIn(string fname, lli x1, lli x2, lli x3, lli x4, lli x5) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << "," << x4 << "," << x5 <<")" << endl;}

void debugOut(string fname, lli x1, lli x2, lli x3, lli ans) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << ") = " << ans << endl;}

void debugOut(string fname, lli x1, lli x2, lli x3, lli x4, lli ans) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << "," << x4 << ") = " << ans << endl;}

void debugOut(string fname, lli x1, lli x2, lli x3, lli x4, lli x5, lli ans) {
	cout << string(rec_depth, ' ') << fname << "(" << x1 << "," << x2 <<	"," << x3 << "," << x4 << "," << x5 << ") = " << ans << endl;}



// fill table for binomial coefficients.
void chooseInit(int n) {
	table_choose.resize(n+1, vector<lli> (n+1, 0));
	table_choose[0][0] = 1;
	
	for(int m = 1; m <= n; ++m) {
		table_choose[m][0] = 1;
		for(int k = 1; k <= m; ++k) {
			table_choose[m][k] = table_choose[m-1][k] + table_choose[m-1][k-1];
		}
	}
}

// 10 choose 3 = 10*9*8 / 1 * 2 * 3
// n choose k for n > k defaults to 0.
lli choose(int n, int k) {
	if(n < k) return 0;
	return table_choose[n][k];
}

// Counts connected chordal graphs on k vertices
// assume k > 0
lli chordal_conn(int k) {
	rec_depth++;
	
	if (table_chordal_conn[k] == notSaved) table_chordal_conn[k] = comp_chordal_conn(k);
	lli ans = table_chordal_conn[k];
	
	rec_depth--;

	return ans;
}

lli comp_chordal_conn(int k) {
	//empty graph counts here
	// clique evap at time 1 countrs here for t>=2, not for t=1
	lli ans = 0;
		
	// clique evap at time 1 countrs here for t=1, not higher t. 
	for(int t = 1; t <= k; ++t) {
		for(int l = 1; l <= k; ++l) {
			ans += choose(k,l) * f(t, 0, l, k-l);
		}
	}

	return ans;
}

// Counts chordal graphs on k vertices
// assume k > 0
lli chordal(int k) {		
	rec_depth++;
	
	if (table_chordal[k] == notSaved) table_chordal[k] = comp_chordal(k);
	lli ans = table_chordal[k];
	
	rec_depth--;

	return ans;
}

lli comp_chordal(int k) {
	// Base case k = 0
	if(k == 0) return 1;

	// Guess the set of vertices in the first component
	lli ans = 0;
	for (int kk = 1; kk <= k; ++kk) {
		ans += choose(k - 1, kk - 1) * chordal_conn(kk) * chordal(k - kk);
	}
	
	return ans;
}

// chordal graphs that evaporate at or before t
// x special vertices (do NOT allow 0)
// k normal vertices
lli g(int t, int x, int z, int k) {
	
	if(verbatim) debugIn("g",t,x,z,k);
		
	rec_depth++;
	
	if(table_g[t][x][z][k] == notSaved) table_g[t][x][z][k] = comp_g(t, x, z, k);
	lli ans = table_g[t][x][z][k];
	
	rec_depth--;
	
	if(verbatim) debugOut("g",t,x,z,k,ans);
	
	return ans;

}

lli comp_g(int t, int x, int z, int k) {
	// Base case t=0
	if(t == 0) return k == 0 ? 1 : 0;
		
	// non zero t. If x is zero cant use x for connectivity, so have
	// to ensure connectivity some other way. Call f maybe??
	if(x == 0) {
		// EXCEPTION HERE
		exit(1);
	}	
		
	// main case: t > 0 and x > 0
	// Guess the set of vertices in components evaporate at time
	// Precisely t.	
	lli ans = 0;
	for(int kk = 0; kk <= k; ++kk) {
		ans += choose(k, kk) * gTilde(t, x, z, kk) * g(t - 1, x, z, k - kk);
	}
	
	return ans;
}

// Same as g, but now all connected components evaporate at time exactly t
// graphs where V = X count (vacously true)
//
// Disallowing x=0
//
// Guess component that contains lowest labeled vertex.
lli gTilde(int t, int x, int z, int k) {
	
	if(verbatim) debugIn("gTilde", t, x, z, k);
	
	rec_depth++;
	
	if(table_gTilde[t][x][z][k] == notSaved) table_gTilde[t][x][z][k] = comp_gTilde(t, x, z, k);
	lli ans = table_gTilde[t][x][z][k];
	
	rec_depth--;
	
	if(verbatim) debugOut("gTilde",t,x,z,k,ans);
	
	return ans;
	
}

lli comp_gTilde(int t, int x, int z, int k) {
	// base case: k=0
	if(k == 0) return 1;
	
	// if x=0 we cant use X for connectivity. Either define this
	// or remove x=0 from domain 
	if(x == 0) {
		// EXCEPTION HERE
		exit(1);
	}
	
	lli ans = 0;
	
	// t and x both nonzero.
	// guess #labels in comp.
	// guess its neigh into x
	for(int kk = 1; kk <= k; ++kk) {
		for(int xx = 1; xx <= x; ++xx) {
			ans += choose(k-1, kk-1) * (choose(x, xx) - choose(z, xx)) * g1Tilde(t, xx, kk) * gTilde(t, x, z, k-kk);
		}
	}
	
	return ans;
	
}

// Same as gTilde, but now no component of G-X sees all of X
// graphs where V = X count (vacously true)
//
// Disallowing x=0
//
// Guess component that contains lowest labeled vertex.
lli gHat(int t, int x, int z, int k) {
	
	if(verbatim) debugIn("gHat",t,x,z,k);
	
	rec_depth++;
	
	if(table_gHat[t][x][z][k] == notSaved) table_gHat[t][x][z][k] = comp_gHat(t, x, z, k);
	lli ans = table_gHat[t][x][z][k];

	rec_depth--;
	
	if(verbatim) debugOut("gHat",t,x,z,k,ans);
	
	return ans;
}

lli comp_gHat(int t, int x, int z, int k) {

	// base case: k=0
	if(k == 0) return 1;
	
	// if x=0 we cant use X for connectivity. Either define this
	// or remove x=0 from domain 
	if(x == 0) {
		// EXCEPTION HERE
		exit(1);
	}
	
	lli ans = 0;
	
	// t and x both nonzero.
	// guess #labels in comp.
	// guess its neigh into x
	// now xx < x because no component sees all of x.
	for(int kk = 1; kk <= k; ++kk) {
		for(int xx = 1; xx < x; ++xx) {
			ans += choose(k-1, kk-1) * (choose(x, xx) - choose(z, xx)) * g1Tilde(t, xx, kk) * gHat(t, x, z, k-kk);
		}
	}
	
	return ans;
}

// exactly one comp of G-X, evap at exactly t. Sees all of X.
// Disallowing x=0 
lli g1Tilde(int t, int x, int k) {
	
	if(verbatim) debugIn("g1Tilde", t, x, k);
	
	rec_depth++;
	
	if(table_g1Tilde[t][x][k] == notSaved) table_g1Tilde[t][x][k] = comp_g1Tilde(t, x, k);
	lli ans = table_g1Tilde[t][x][k];
	
	rec_depth--;
	
	if(verbatim) debugOut("g1Tilde", t, x, k, ans);
	
	return ans;
}

lli comp_g1Tilde(int t, int x, int k) {
	//base case k=0 return 0
	if(k == 0) return 0;
	
	// nonzero k, t=0, also 0
	if(t == 0) return 0;
	
	// nonzero t and k, x zero. Shoulf be undefined?? 
	if(x == 0) {
		// EXCEPTION
		exit(1);
	}
	
	lli ans = 0;
	
	// guess the label set of L
	for(int l = 1; l <= k; ++l) {
		ans += choose(k, l) * f(t, x, l, k - l);
	}
	
	return ans;
}

// at least two comp of G-X, evap at exactly t. Sees all of X.
// Disallowing x=0 
lli g2Tilde(int t, int x, int k) {
		
		if(verbatim) debugIn("g2Tilde", t, x, k);
		
		rec_depth++;
		
		if(table_g2Tilde[t][x][k] == notSaved) table_g2Tilde[t][x][k] = comp_g2Tilde(t, x, k);
		lli ans = table_g2Tilde[t][x][k];
		
		rec_depth--;
		
		if(verbatim) debugOut("g2Tilde", t, x, k, ans);
		
		return ans;
}

lli comp_g2Tilde(int t, int x, int k) {
	//base case k=0 return 0
	if(k == 0) return 0;
	
	// nonzero k, t=0, also 0)
	if(t == 0) return 0;
	
	// nonzero t and k, x zero. Shoulf be undefined?? 
	if(x == 0) {
		// EXCEPTION
		exit(1);
	}
	
	lli ans = 0;
	
	// guess the label set of comp of G-X with lowest label
	for(int kk = 1; kk < k; ++kk) {
		ans += choose(k - 1, kk - 1) * g1Tilde(t, x, kk) * ( g1Tilde(t, x, k-kk) + g2Tilde(t, x, k-kk) );
	}
	
	return ans;
}

// Connected chordal graphs evap at exactly t, |X|=x, L=l, |rest|=k
// Require G-X connected and L union X a clique. 
// allow x = 0
// empty graph (or empty nonX -- i.e l+k=0) evaps at time t=0, otherwise l must be non zero. If t==1 then k must be 0, for t>= 2 k must be non zero.
lli f(int t, int x, int l, int k) {
		
		if(verbatim) debugIn("f", t, x, l, k);
		
		rec_depth++;
		
		if(table_f[t][x][l][k] == notSaved) table_f[t][x][l][k] = comp_f(t, x, l, k);
		lli ans = table_f[t][x][l][k];
		
		rec_depth--;
		
		if(verbatim) debugOut("f", t, x, l, k, ans);
		
		return ans;
}

lli comp_f(int t, int x, int l, int k) {
	// If size of clique is too large, no dice.
	if(x + l > w) return 0;

	// base case, t=0
	if(t == 0) return 0;
	// if(t == 0) return l+k==0 ? 1 : 0;
	
	// additional base case, t=1. Then L must be non zero and k must be zero. x can be anything, including 0 
	if(t == 1) return k == 0 ? 1 : 0;
	
	// t >= 2 then k must be non zero.
	if(k == 0) return 0;
	
	// t >= 2, k >= 1, l >= 1
	lli ans = 0;
	
	// guess the set of labels in components that evaporate at time exactly t-1. 
	// this set must be non zero.
	for(int kk = 1; kk<=k; ++kk) {
		ans += choose(k, kk) * fTilde(t, x, l, kk) * g(t-2, x+l, x, k-kk);
	}
	
	return ans;
}

// Same as f but now every component of G-(X u L) evaporates at time exactly t-1
// Require at least one such component!
// Require G-X connected and L union X a clique. 
// allow x = 0
lli fTilde(int t, int x, int l, int k) {

	if(verbatim) debugIn("fTilde", t, x, l, k);

	rec_depth++;
	
	if(table_fTilde[t][x][l][k] == notSaved) table_fTilde[t][x][l][k] = comp_fTilde(t, x, l, k);
	lli ans = table_fTilde[t][x][l][k];
	
	rec_depth--;

	if(verbatim) debugOut("fTilde", t, x, l, k, ans);

	return ans;
}

lli comp_fTilde(int t, int x, int l, int k) {
	// base case, t=0.
	if(t == 0) return 0;
	
	// second base case t=1, then comps of G-L union X evap at time 0, this is not possible.
	if(t == 1) return 0;
	
	// since t>= 2, k must be non zero.
	if(k == 0) return 0;	
	
	// t >= 2, k >= 1, l >= 1
	lli ans = 0;
	
	ans += fHat(t, x, l, k);
	
	// precisely one comp of G-(X u L) sees all of X u L
	// since l is non zero, at least one vertex is not in the comp, so kk < k
	for(int kk = 1; kk < k; ++kk) {
		ans += choose(k, kk) * g1Tilde(t-1, x+l, kk) * fHat(t,x,l,k-kk);
	}

	// at least 2 comps of G-(X u L) see all of X u L.
	// at least 1 vertex in those componnets. Possibly all vertices. 
	for(int kk = 1; kk <= k; ++kk) {
		ans += choose(k, kk) * g2Tilde(t - 1, x + l, kk) * gHat(t - 1, x + l, x, k - kk);
	}

	return ans;
}

// Same as fTilde but now no component sees all of X union L
// Require at least one such component of G - (X union L)
// Require G-X connected and L union X a clique. 
// allow x = 0
lli fHat(int t, int x, int l, int k) {
	
	if(verbatim) debugIn("fHat", t, x, l, k);
	
	rec_depth++;
	
	if(table_fHat[t][x][l][k] == notSaved) table_fHat[t][x][l][k] = comp_fHat(t, x, l, k);
	lli ans = table_fHat[t][x][l][k];

	rec_depth--;
	
	if(verbatim) debugOut("fHat", t, x, l, k, ans);
	
	return ans;
}

lli comp_fHat(int t, int x, int l, int k) {
	return fHat(t, x, x, l, k);
}

// Overload of fHat but now with z -- every comp of G-X must see someone outside of 1...z
lli fHat(int t, int x, int z, int l, int k) {
	
	if(verbatim) debugIn("fHat", t, x, z, l, k);
	
	rec_depth++;
	
	if(table_fHat5[t][x][z][l][k] == notSaved) table_fHat5[t][x][z][l][k] = comp_fHat(t, x, z, l, k);
	lli ans =  table_fHat5[t][x][z][l][k];
	
	rec_depth--;
	
	if(verbatim) debugOut("fHat", t, x, z, l, k, ans);
	
	return ans;
}

lli comp_fHat(int t, int x, int z, int l, int k) {
	
	// base case, t=0.
	if(t == 0) return 0;
	
	// second base case t=1, then comps of G-L union X evap at time 0, this is not possible. 
	if(t == 1) return 0;
	
	// since t>= 2, k must be non zero.
	if(k == 0) return 0;
	
	// t >= 2, k >= 1, l >= 1
	lli ans = 0;
	
	// guess component comtaining the lowest labeled vertex
	for(int kk = 1; kk <= k; ++kk) {
		for(int xx = 0; xx <= x; ++xx) {
			for(int ll = 0; ll <= l; ++ll) {
				if(xx+ll == 0 || xx+ll == x+l) continue;
				
				lli prod = 1;
				prod *= choose(k-1, kk-1) * choose(l, ll);
				// prod *= (ll > 0 || z == 0) ? choose(x, xx) : choose(x, xx) - choose(z, xx);
				prod *= (ll > 0) ? choose(x, xx) : choose(x, xx) - choose(z, xx);
				prod *= g1Tilde(t-1, xx + ll, kk);
				prod *= ll < l ? fHat(t, x+ll, z, l-ll, k-kk) : gHat(t-1, x+ll, z, k-kk );
				
				ans += prod;
				
			}
		}
	}

	return ans;

}

int main() {
	int n;
	cout << "First input is the number of vertices, second is the bound on the clique size." << endl;
	cin >> n >> w;
	
	chooseInit(n);

	table_chordal.resize(n+1, notSaved);
	table_chordal_conn.resize(n+1, notSaved);
	table_g.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_gTilde.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_gHat.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_g1Tilde.resize(n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved)));
	table_g2Tilde.resize(n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved)));
	table_f.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_fTilde.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_fHat.resize(n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved))));
	table_fHat5.resize(n+1, vector<vector<vector<vector <lli> > > > (n+1, vector<vector<vector<lli> > > (n+1, vector<vector<lli> > (n+1, vector<lli> (n+1, notSaved)))));
	
	cout << "number of connected labeled chordal graphs: " << chordal_conn(n) << endl;
	cout << "number of labeled chordal graphs: " << chordal(n) << endl;
}
