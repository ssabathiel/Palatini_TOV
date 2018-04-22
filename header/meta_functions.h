#ifndef META_FUNCTIONS
#define META_FUNCTIONS

typedef std::pair<double, double> (*get_gradients_functions) (double m, double press, double r);
typedef double (*p_of_rho_functions) (double rho);
typedef double (*rho_of_p_functions) (double p);

typedef double (*drho_dp_functions) (double p);
typedef double (*ddrho_dpp_functions) (double p);

extern std::vector<p_of_rho_functions> rho_of_p;


#endif // META_FUNCTIONS

