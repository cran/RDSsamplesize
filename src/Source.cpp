//[[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>
#include <functional>

using namespace std;
using namespace Rcpp;

struct At {
    // structure to collect A_t results
    int zk;
    int nk;
    double g;
};
struct Ft {
    // structure to collect F_t results
    int zk;
    double pmf;
};



bool IsNAN(double i) { return (std::isnan(i)); }

double binom_coeff(int n, int k) {
    //Binomial coefficient evaluation.
    return(exp(lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1)));
}
void dist_w1(vector<struct At> &A1, vector<struct Ft> &F1, int n0, int m, double p) {
   //Get distribution of the number of recruits by the 1st wave (Z1).
    for (int n1 = 0; n1 <= (n0 * m);n1++) {
        At temp;
        temp.zk = n1;
        temp.nk = n1;
        temp.g = binom_coeff(n0 * m, n1) * pow(p,n1)  * pow(1 - p,n0 * m - n1);
        A1.push_back(temp);
    }
    F1.resize(n0 * m + 1);
    for (const auto& s : A1) {
        F1[s.zk].zk = s.zk;
        F1[s.zk].pmf += s.g;
    }
    return;
}

double func(double n0, int n1, int m, double p) {
    //Stirling's approximation of the binomial function
    return (n0*m+.5)*log(n0*m)-(n0*m-n1+.5)*log(n0*m-n1)+n0*m*log(1-p);
}

double Df(double n0, int n1, int m, double p) {
    //Derivative of the "func" above
    return m * log(n0 * m) + (n0 * m + 0.5) * (m/(n0 * m)) - (m * log(n0 *
    m - n1) + (n0 * m - n1 + 0.5) * (m/(n0 * m - n1))) + m *
    log(1 - p);
}

int bisection(int n1, int m, double p, int L, int U) {
   // Bisection method to find the maximizer given a region.
   double Df_low;
   do {
       Df_low = Df(L, n1, m, p);
       L++;
   } while (IsNAN(Df_low));
   double Df_up = Df(U, n1, m, p);
   if ( Df_low* Df_up >= 0) {
      double f_low = func(L, n1, m, p);
      double f_up = func(U, n1, m, p);
      if(f_low<f_up)
          return U;
      else
          return L;
   }
   double c = L;
   while (U-L > 1) {
      c = (U+L)/2;
      if (Df(c, n1, m, p) == 0.0)
         break;
      else if (Df(c, n1, m, p)*Df(L, n1, m, p) < 0)
         U = c;
      else
         L = c;
   }
   return floor(c);
}

int binom_max(int n1, int m, double p, int L, int U) {
    //find the maximizer for a binomial function given a region.
    if (L == U)
        return U;
    int maximizer = bisection(n1, m, p, L, U);
    return maximizer;
}
double binom_eval(int n1, int n0, int m, double p) {
    // Binomial function evaluation.
    return exp(lgamma(n0 * m + 1) - lgamma(n1 + 1) - lgamma(n0 * m - n1 + 1)+n1*log(p)+
    (n0*m-n1)*log(1-p));
}

double gk(int k, int m, double p, int nk, int n, double thres, vector<At>prev){
    //Calculation of gk(nk,n):=Pr(nk new recruits at k-th wave & a total of n recruits by k-th wave & ongoing recruitment process).

    auto it = max_element(prev.begin(), prev.end(),
        [](const At& a, const At& b) { return a.nk < b.nk; });
    if (it == prev.end()) throw "max_element called on emtpy vector";
    auto max_nk = it->nk;
    if (nk > max_nk * m)
        return 0.0;

    double out = 0.0;
    vector<int> l1;
    vector<At>::iterator it_l0 = prev.begin();
    while ((it_l0 = find_if(it_l0, prev.end(), [&](const At& s) {return s.g > thres; })) != prev.end())
    {
        int loc=distance(prev.begin(), it_l0);
        if (prev[loc].nk >= 1 + ((nk - 1) / m))
            l1.push_back(prev[loc].nk);
        it_l0++;
    }
    if (l1.size() == 0)
        return 0.0;
    int lower=*min_element(l1.begin(), l1.end());
    int upper = *max_element(l1.begin(), l1.end());
    int maximizer = binom_max(nk, m, p, lower, upper);

    int left = maximizer;
    double b = binom_eval(nk, left, m, p);

    while ((left >= lower) & (b > thres)) {
        auto findnot = find(l1.begin(), l1.end(), left);
        if (findnot != l1.end()) {
            auto found = find_if(prev.begin(), prev.end(),
                [&](const At& s) {return s.nk == left; });
            out += b * found->g;
        }
        left += -1;
        b = binom_eval(nk, left, m, p);
    };
    int right = maximizer + 1;
    if (right > upper)
        return out;
    b = binom_eval(nk, right, m, p);
    while ((right <= upper) & (b > thres)) {
        auto findnot = find(l1.begin(), l1.end(), right);
        if (findnot != l1.end()) {
            auto found = find_if(prev.begin(), prev.end(),
                [&](const At& s) {return s.nk == right; });
            out += b * found->g;
        }
        right += 1;
        b = binom_eval(nk, right, m, p);
    };
    return out;
}

namespace Utils {
    template<typename InputIt, typename T = typename std::iterator_traits<InputIt>::value_type, typename Pred>
    std::vector<T> distinct(InputIt first, InputIt last, const Pred& predicate) {
        std::vector<T> results;
        std::copy_if(first, last, std::back_inserter(results),
            [&](const T& lhs) {
                return std::find_if(results.begin(), results.end(), [&](const T& rhs) {
                    return predicate(lhs, rhs);
                    }) == results.end();
            }
        );
        return results;
    }
    template<typename InputIt, typename T = typename std::iterator_traits<InputIt>::value_type>
    std::vector<T> distinct(InputIt first, InputIt last) {
        auto pred = [](const T& x, const T& y) { return x == y; };
        return distinct(first, last, pred);
    }
}

void gk_sum_nk(vector<At>&A_n,
    int k, int m, double p, int n, double thres, vector<At>prev) {
   //Calculate sum of gk(nk,n) over nk, equivalently, Pr(n recruits by k-th wave & ongoing recruitment process), and keep documentation

    auto it = max_element(prev.begin(), prev.end(),
        [](const At& a, const At& b) { return a.nk < b.nk; });
    if (it == prev.end()) throw "max_element called on emtpy vector";
    auto max_nk = it->nk;
    it = max_element(prev.begin(),prev.end(),
        [](const At& a, const At& b) { return a.zk < b.zk; });
    if (it == prev.end()) throw "max_element called on emtpy vector";
    auto max_zk = it->zk;
    int max = max_nk * m + max_zk;

    if ((n > max) | (n < k)) {
        return;
    }
    auto find_uniq_zk = [](const At& lhs, const At& rhs){return lhs.zk == rhs.zk; };
    auto results = Utils::distinct(prev.begin(), prev.end(), find_uniq_zk);
    for (const auto& uniq : results) {
        int nk = n - uniq.zk;
        if (nk <= 0)
            continue;
        vector<At> sub;
        vector<At>::iterator it = prev.begin();
        while ((it = find_if(it, prev.end(), [&](const At& s) {return s.zk == uniq.zk; })) != prev.end())
        {
            int loc = distance(prev.begin(), it);
            sub.push_back(prev[loc]);
            it++;
        };

        At temp;
        temp.g = gk(k, m, p, nk, n, thres, sub);
        temp.zk = n;
        temp.nk = nk;
        if(temp.g>0)
            A_n.push_back(temp);
    }
    if (thres > accumulate(A_n.begin(), A_n.end(), 0.0, [](double sum, const At& elem) { return sum + elem.g; })) {
        A_n.clear();
        return;
    }

    return;
}

void pk(vector<struct At> &A_n,double &pk_alive,
    int k, int m, vector<double>p_list, int n, vector<double>thres_list,
    vector<vector<At>>Aall) {
 //  Calculation of Pr(n recruits by k-th wave) and keep documentation of exact cases for ongoing recruitment process.

    gk_sum_nk(A_n,k, m, p_list[k - 1], n, thres_list[k - 2], Aall[k - 2]);
    pk_alive = accumulate(A_n.begin(), A_n.end(), 0.0, [](double sum, const At& elem) { return sum + elem.g; });

    int t = k;
    while (t > 1) {
        auto it = max_element(Aall[t - 2].begin(), Aall[t - 2].end(),
            [](const At& a, const At& b) { return a.zk < b.zk; });
        if (it == Aall[t - 2].end()) throw "max_element called on emtpy vector";
        auto max_zk = it->zk;
        if (n > max_zk)
            break;

        vector<At>::iterator it_loc = Aall[t - 2].begin();
        while ((it_loc = find_if(it_loc, Aall[t - 2].end(), [&](const At& s) {return s.zk == n; })) != Aall[t - 2].end())
        {
            int loc = distance(Aall[t - 2].begin(), it_loc);
            pk_alive += Aall[t - 2][loc].g * pow((1 - p_list[t - 1]), Aall[t - 2][loc].nk * m);
            it_loc++;
        }
        t--;
    }

    return;
}

void dist_wk(vector<struct At>&Ak, vector<struct Ft>&Fk,
    vector<vector<struct At>>Aall, vector<struct Ft> F_prev,
    int k, int m, vector<double> p_list, vector<double>thres_list) {
    //Get distribution of the number of recruits by the k-th wave (Zk).

    vector<At> A_prev = Aall.back();
    auto it = max_element(A_prev.begin(), A_prev.end(),
        [](const At& a, const At& b) { return a.nk < b.nk; });
    if (it == A_prev.end()) throw "max_element called on emtpy vector";
    auto max_nk = it->nk;

    it = max_element(A_prev.begin(), A_prev.end(),
        [](const At& a, const At& b) { return a.zk < b.zk; });
    if (it == A_prev.end()) throw "max_element called on emtpy vector";
    auto max_zk = it->zk;
    int max = max_nk * m + max_zk;

    Ft temp;
    double condSum = 0.0;
    for (int n = 0; n <= max; n++) {
        double pk_alive;
        vector<At> A_n;

       pk(A_n, pk_alive, k, m, p_list, n, thres_list, Aall);
       if ((pk_alive == 0) & (condSum > .5))
            break;
       if (A_n.size() > 0)
           Ak.insert(Ak.end(),A_n.begin(),A_n.end());
        temp.zk = n;
        temp.pmf = pk_alive;
        condSum += pk_alive;
        if (pk_alive > 0)
            Fk.push_back(temp);
    }
    return;
}

double q(int m, double p, double s) {
    return pow(p * s + 1 - p, m);
}

double G_nt(int m, int n0, vector<double>p_list, double s, int t) {
    //PGF calculation
    if (t == 0)
        return pow(s, n0);
    else//if (t > 0)
        return G_nt(m, n0, p_list, q(m, p_list[t - 1], s), t - 1);
}

vector<double>P_tau(int m, int n0, vector<double>p_list, int maxT) {
    //Extinction probability calculation
    double s = 0;
    vector<double>G_nt_list;
    for (int t = 0; t <= maxT; t++) {
        G_nt_list.push_back(G_nt(m, n0, p_list, s, t));
    }
    adjacent_difference(G_nt_list.begin(), G_nt_list.end(), G_nt_list.begin());
    G_nt_list.erase(G_nt_list.begin());
    return G_nt_list;
}
int check_maxT(vector<double>P_tau_list, double cutoff) {
    //  Check the validity of user's input maxT.
    vector<double>cumsum(P_tau_list.size());
    partial_sum(P_tau_list.begin(), P_tau_list.end(), cumsum.begin(), plus<double>());
    auto it_maxT = find_if(cumsum.begin(), cumsum.end(),
        [&](const double& val) { return val >= cutoff; });
    if (it_maxT != cumsum.end())
        return it_maxT - cumsum.begin()+1;
    else
        return P_tau_list.size();
}
vector<double>ind_tol(vector<double>P_tau_list, int maxT_checked, double tol) {
    //  Calculate the tolenrance value of accuracy loss at each wave.
    vector<double>P_tau_cut(P_tau_list.begin(), P_tau_list.begin() + maxT_checked);
    double sum = accumulate(P_tau_cut.begin(), P_tau_cut.end(), 0.0);
    transform(P_tau_cut.begin(), P_tau_cut.end(),P_tau_cut.begin(), [&](double& p) {return p * tol / sum; });
    return P_tau_cut;
}
double thres(vector<double>x, double ind_tol_val) {
    //Threshold calculation
    sort(x.begin(), x.end());

    vector<double>cumsum(x.size());
    partial_sum(x.begin(), x.end(), cumsum.begin(), plus<double>());
    auto it_maxT = find_if(cumsum.begin(), cumsum.end(),
        [&](const double& val) { return val > ind_tol_val; });
    int index = 0;
    if (it_maxT != cumsum.end()) {
        index = it_maxT - cumsum.begin();
        if (index - 1 > 0)
            return x[index - 1];
        else
            return x[0];
    }
    else
        return x[index];
}



void sample_size_eval(vector<vector<At>>& Aall, vector<vector<Ft>>& Fall,
    int n0, int m, int maxT, vector<double>p_list, double tol) {
    //Sample size power calculation for each wave.
    vector<double>thres_list;
    vector<double>P_tau_list = P_tau(m, n0, p_list, maxT);
    double cutoff = .95;
    maxT = check_maxT(P_tau_list, cutoff);
    vector<double>ind_tol_list = ind_tol(P_tau_list, maxT, tol);

    //1st wave
    vector<struct At>A1;
    vector<struct Ft>F1;
    dist_w1(A1, F1, n0, m, p_list[0]);
    Aall.push_back(A1);
    Fall.push_back(F1);


    //after 1st wave
    for (int k = 2; k <= maxT; k++) {
        vector<double> x(Aall.back().size());
        transform(Aall.back().begin(), Aall.back().end(), x.begin(),
            mem_fn(&At::g));
        double ind_thres = thres(x, ind_tol_list[k - 2]);
        thres_list.push_back(ind_thres);

        vector<struct At>Ak;
        vector<struct Ft>Fk;
        dist_wk(Ak, Fk, Aall, Fall.back(), k, m, p_list, thres_list);

        Aall.push_back(Ak);
        Fall.push_back(Fk);
        if ((Ak.size() == 0) | (Fk.size() == 0)){
            Rcpp::Rcout<<"break ahead "<<Ak.size()<<" "<<Fk.size() ;
            break;
        }
        Ak.clear();
        Fk.clear();
    }

    return;
}

//[[Rcpp::export]]
List size(int n0, int m, int maxT, NumericVector p_list_vec, double tol){
    //main function
    //read in p_list
    vector<double>p_list = as<vector<double>> (p_list_vec);
    vector<vector<At>> Aall;
    vector<vector<Ft>> Fall;
    sample_size_eval(Aall, Fall,n0, m, maxT,p_list, tol);

    List output(Fall.size()+1);
    vector<double>P_tau_list = P_tau(m, n0, p_list, maxT);
    output(0)=P_tau_list;

    for (int i = 0; i < Fall.size(); i++) {
      vector<Ft> F_n=Fall[i];
      NumericMatrix F_temp(F_n.size(),2);
      for (int j = 0; j < F_n.size(); j++) {
        F_temp(j,0)=F_n[j].zk;
        F_temp(j,1)=F_n[j].pmf;
      }
      output(i+1)=F_temp;

    }

    return output;
}
