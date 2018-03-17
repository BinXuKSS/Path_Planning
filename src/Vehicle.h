#ifndef VEHICLE
#define VEHICLE

#include <vector>
#include <map>
#include <string>

using namespace std;

Class Vehicle
{
	public:
		double s;
		double s_d;
		double s_dd;
		double d;
		double d_d;
		double d_dd;

		string state;

		vector<string> possible_state;
		vector<double> s_trajectory_coeffs, d_trajectory_coeffs;

		/* Constructors */
		Vehicle();
		Vehicle(double s, double s_d, double s_dd, double d, double d_d, double d_dd);

		/* Destructor */
		
		virtual ~Vehicle();
		
		
}
#endif


