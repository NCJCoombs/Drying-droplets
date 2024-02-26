//====================================================================
// Driver for a 1D drying drop problem
//====================================================================

#include <fenv.h>

// Generic routines
#include "generic.h"

// The equations
#include "thin_film_brinkman.h"

// The mesh
#include "meshes/one_d_mesh.h"

using namespace std;
using namespace oomph;


/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////

//======start_of_GlobalPhysicalVariables==============================
/// Namespace for global physical variables
//====================================================================
namespace GlobalPhysicalVariables
{
  /// Capillary number
  double Ca = 1.0e-5;

  /// Jamming threshold (random close packing)
  double phi_c = 0.64;

  // Initial solute concentration
  double phi_0 = 0.0235; // phi_c / 25.0;

  /// Peclet number
  double Pe = 200.0;

  /// Scaled inverse pore size
  double nu = 1.0e5;

  /// Axisymmetric geometry?
  bool axisym_flag = true;

  /// Diffusive evaporation?
  unsigned is_diffusive = 0;

  /// The dry-out time: equal to 0.5 regardless of the evaporation mode
  double t_f = 0.5;

  // Parameters for the transition function
  double phi_tilde = 0.7 * phi_c;
  double delta = 0.001;
  double alpha;
  double beta;

  /// Evaporation flux function
  void get_evaporation_flux_fct(const double& time,
                                const Vector<double>& x,
                                double& evap)
  {
    if (is_diffusive)
    {
      evap = 0.5 * pow(1.0 - x[0] * x[0], -0.5);
    }
    else
    {
      evap = 1.0;
    }
  }

  /// Transition function for particle jamming
  void get_transition_fct(const double& phi, double& H)
  {
    alpha = 1.0 -
            2.0 * pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 2.0) +
            pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 4.0);

    beta = 4.0 *
           ((phi_c - phi_tilde) / (phi_c + delta - phi_tilde) -
            pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 3.0)) /
           (alpha * (phi_c + delta - phi_tilde));

    if (phi < phi_tilde)
    {
      H = 1.0;
    }
    else if (phi > phi_tilde && phi < phi_c)
    {
      H = 1.0 -
          2.0 * pow((phi - phi_tilde) / (phi_c + delta - phi_tilde), 2.0) +
          pow((phi - phi_tilde) / (phi_c + delta - phi_tilde), 4.0);
    }
    else
    {
      H = alpha * exp(beta * (phi_c - phi));
    }
  }

  /// Derivative of transition function for particle jamming
  void get_transition_fct_deriv(const double& phi, double& dHdphi)
  {
    alpha = 1.0 -
            2.0 * pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 2.0) +
            pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 4.0);

    beta = 4.0 *
           ((phi_c - phi_tilde) / (phi_c + delta - phi_tilde) -
            pow((phi_c - phi_tilde) / (phi_c + delta - phi_tilde), 3.0)) /
           (alpha * (phi_c + delta - phi_tilde));

    if (phi < phi_tilde)
    {
      dHdphi = 0.0;
    }
    else if (phi > phi_tilde && phi < phi_c)
    {
      dHdphi = -4.0 *
               ((phi - phi_tilde) / (phi_c + delta - phi_tilde) -
                pow((phi - phi_tilde) / (phi_c + delta - phi_tilde), 3.0)) /
               (phi_c + delta - phi_tilde);
    }
    else
    {
      dHdphi = -alpha * beta * exp(beta * (phi_c - phi));
    }
  }

  /// Mobility term
  void get_mobility(const double& h, const double& phi, double& Q)
  {
    double mu_KD_inv;
    if (phi < phi_c)
    {
      mu_KD_inv = pow(1.0 - phi / phi_c, 1.6);
    }
    else
    {
      mu_KD_inv = 0.0;
    }

    Q = -mu_KD_inv * pow(h, 3.0) / (3.0 * Ca) - h / (Ca * pow(nu, 2.0));
  }

  /// Derivative of mobility w.r.t. h
  void get_partial_mobility_partial_h(const double& h,
                                      const double& phi,
                                      double& dQdh)
  {
    double mu_KD_inv;
    if (phi < phi_c)
    {
      mu_KD_inv = pow(1.0 - phi / phi_c, 1.6);
    }
    else
    {
      mu_KD_inv = 0.0;
    }

    dQdh = -pow(h, 2.0) * mu_KD_inv / Ca - 1.0 / (pow(nu, 2.0) * Ca);
  }

  /// Derivative of mobility w.r.t. c
  void get_partial_mobility_partial_phi(const double& h,
                                        const double& phi,
                                        double& dQdphi)
  {
    double mu_KD_inv_deriv;
    if (phi < phi_c)
    {
      mu_KD_inv_deriv = -1.6 * pow(1.0 - phi / phi_c, 0.6) / phi_c;
    }
    else
    {
      mu_KD_inv_deriv = 0.0;
    }

    dQdphi = -pow(h, 3.0) * mu_KD_inv_deriv / (3 * Ca);
  }

  void get_solute_mobility(const double& h, const double& phi, double& Q)
  {
    double mu_KD_inv;
    if (phi < phi_c)
    {
      mu_KD_inv = pow(1.0 - phi / phi_c, 1.6);
    }
    else
    {
      mu_KD_inv = 0.0;
    }

    double transition_fct;
    get_transition_fct(phi, transition_fct);

    Q = -transition_fct * mu_KD_inv * pow(h, 3.0) / (3.0 * Ca);
  }

  /// Derivative of solute mobility w.r.t. h
  void get_partial_solute_mobility_partial_h(const double& h,
                                             const double& phi,
                                             double& dQdh)
  {
    double mu_KD_inv;
    if (phi < phi_c)
    {
      mu_KD_inv = pow(1.0 - phi / phi_c, 1.6);
    }
    else
    {
      mu_KD_inv = 0.0;
    }

    double transition_fct;
    get_transition_fct(phi, transition_fct);


    dQdh = -transition_fct * mu_KD_inv * pow(h, 2.0) / Ca;
  }

  /// Derivative of solute mobility w.r.t. phi
  void get_partial_solute_mobility_partial_phi(const double& h,
                                               const double& phi,
                                               double& dQdphi)
  {
    double mu_KD_inv;
    if (phi < phi_c)
    {
      mu_KD_inv = pow(1.0 - phi / phi_c, 1.6);
    }
    else
    {
      mu_KD_inv = 0.0;
    }

    double mu_KD_inv_deriv;
    if (phi < phi_c)
    {
      mu_KD_inv_deriv = -1.6 * pow(1.0 - phi / phi_c, 0.6) / phi_c;
    }
    else
    {
      mu_KD_inv_deriv = 0.0;
    }

    double transition_fct;
    get_transition_fct(phi, transition_fct);

    double transition_fct_deriv;
    get_transition_fct_deriv(phi, transition_fct_deriv);

    dQdphi =
      -(transition_fct_deriv * mu_KD_inv + transition_fct * mu_KD_inv_deriv) *
      pow(h, 3.0) / (3.0 * Ca);
  }
} // namespace GlobalPhysicalVariables


//======start_of_GlobalSimSettings====================================
/// Namespace for gloabal simulation settings
//====================================================================
namespace GlobalSimSettings
{
  /// Result Folder: Go to into Kinetic or Diffusive directories depending on
  /// the evaporation flag
  string result_folder;

  // Initial number of elements
  unsigned n_element = 100;

  // Maximum permitted error target (elements refined if error is above this)
  double max_permitted_error = 1.0e-7;

  // Minimum permitted error target (elements unrefined if error is below this)
  double min_permitted_error = 1.0e-8;

  // Minimum allowed timestep
  double dt_min = 1.0e-14;

  // Target error for adaptive timestepping
  double epsilon_t = 1.0e-3;

  // How frequently the mesh is rebuilt
  unsigned adapt_number = 3;

  // Max. refinement level: default is five
  unsigned max_refinement_level = 7;

  // Maximum simulation time
  double t_max = 0.6;

  // Computational timestep
  double dt_comp = 1.0e-5;
} // namespace GlobalSimSettings


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class ThinFilmProblem : public virtual Problem
{
public:
  /// Constructor
  ThinFilmProblem();

  /// Destructor
  ~ThinFilmProblem(){};

  /// Update the problem specs after timestep (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep (empty)
  void actions_before_implicit_timestep() {}

  /// Set initial condition (incl previous timesteps) according
  /// to specified function.  Note that his overloads the virtual
  /// function in the Problem base class and is therefore executed
  /// automatically to re-assign the initial conditions during the
  /// spatially adaptive solution at the first timestep.
  void set_initial_condition();

  // Actions after adapt: empty
  void actions_after_adapt() {}

  /// Actions before adapt (empty)
  void actions_before_adapt() {}

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve.
  void actions_before_newton_solve() {}

  /// Doc the solution
  void doc_solution();

  /// Global error norm for adaptive time-stepping
  double global_temporal_error_norm();

private:
  /// Doc info object for labeling output
  DocInfo Doc_info;

  /// Pointers to specific mesh
  RefineableOneDMesh<ELEMENT>* My_mesh_pt;

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
ThinFilmProblem<ELEMENT>::ThinFilmProblem()
{
  // Output directory
  Doc_info.set_directory(GlobalSimSettings::result_folder);

  // Output number
  Doc_info.number() = 0;

  // Allocate the timestepper (this constructs the time object as well)
  // Boolean flag enables adaptive timestepping
  add_time_stepper_pt(new BDF<2>(true));

  // Create the mesh
  My_mesh_pt = new RefineableOneDMesh<ELEMENT>(
    GlobalSimSettings::n_element, 1.0, time_stepper_pt());

  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt = new Z2ErrorEstimator;
  My_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set targets for spatial adaptivity
  My_mesh_pt->max_permitted_error() = GlobalSimSettings::max_permitted_error;
  My_mesh_pt->min_permitted_error() = GlobalSimSettings::min_permitted_error;

  // Max refinement level: large enough so that we never have to worry about it
  My_mesh_pt->max_refinement_level() = GlobalSimSettings::max_refinement_level;

  // Store as the problem's one and only mesh
  Problem::mesh_pt() = My_mesh_pt;

  // Minimum timestep for the problem
  minimum_dt() = GlobalSimSettings::dt_min;

  // Set the boundary conditions for problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned nbound = My_mesh_pt->nboundary();
  for (unsigned ibound = 0; ibound < nbound; ibound++)
  {
    // If we're on the RHS boundary
    if (ibound == 1)
    {
      unsigned num_nod = My_mesh_pt->nboundary_node(ibound);
      for (unsigned inod = 0; inod < num_nod; inod++)
      {
        // Get node
        Node* nod_pt = My_mesh_pt->boundary_node_pt(ibound, inod);

        // If we're on the outer boundary, pin the free surface height only.
        // Leaving the other variable un-pinned imposes a Neumann condition on
        // the flux
        nod_pt->pin(0);
      }
    }
  } // end loop over boundaries


  // Complete the build of all elements so they are fully functional
  for (unsigned e = 0; e < GlobalSimSettings::n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(My_mesh_pt->element_pt(e));

    // Set the evaporation flux pointer
    el_pt->evap_flux_pt() = &GlobalPhysicalVariables::get_evaporation_flux_fct;

    // Set the transition function pointer
    el_pt->transition_fct_pt() = &GlobalPhysicalVariables::get_transition_fct;

    // Set the transition function derivative pointer
    el_pt->transition_fct_deriv_pt() =
      &GlobalPhysicalVariables::get_transition_fct_deriv;

    // Set the mobility function pointer
    el_pt->mobility_fct_pt() = &GlobalPhysicalVariables::get_mobility;

    // Set the mobility derivative dQdh pointer
    el_pt->mobility_fct_deriv_h_pt() =
      &GlobalPhysicalVariables::get_partial_mobility_partial_h;

    // Set the mobility derivative dQdc pointer
    el_pt->mobility_fct_deriv_c_pt() =
      &GlobalPhysicalVariables::get_partial_mobility_partial_phi;

    // Set the solute mobility function pointer
    el_pt->solute_mobility_fct_pt() =
      &GlobalPhysicalVariables::get_solute_mobility;

    // Set the Peclet number pointer
    el_pt->peclet_pt() = &GlobalPhysicalVariables::Pe;

    // Set the dry-out time pointer
    el_pt->t_f_pt() = &GlobalPhysicalVariables::t_f;

    // Set the axisymmetry flag
    el_pt->axisymmetry_flag_pt() = &GlobalPhysicalVariables::axisym_flag;
  }

  // Setup equation numbering scheme
  oomph_info << "Number of equations: " << this->assign_eqn_numbers()
             << std::endl;

} // end_of_constructor


//==start_of_set_initial_condition========================================
/// Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template<class ELEMENT>
void ThinFilmProblem<ELEMENT>::set_initial_condition()
{
  // Set the timestep. Also sets weights on timesteppers
  initialise_dt(GlobalSimSettings::dt_comp);

  // Assign initial condition and do an impulsive start
  // Find number of nodes in mesh
  unsigned num_nod = My_mesh_pt->nnode();

  Vector<double> x(1);

  Node* nod_pt;

  // Loop over the nodes to set initial guess everywhere
  for (unsigned jnod = 0; jnod < num_nod; jnod++)
  {
    // Get node
    nod_pt = My_mesh_pt->node_pt(jnod);

    x[0] = nod_pt->x(0);

    // Dome shape for a perfect circlar domain
    mesh_pt()->node_pt(jnod)->set_value(0, 1.0 - x[0] * x[0]);

    // Initial solute concentration is 1
    mesh_pt()->node_pt(jnod)->set_value(2, GlobalPhysicalVariables::phi_0);
  }

  assign_initial_values_impulsive();
} // End of set_initial_condition


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ThinFilmProblem<ELEMENT>::doc_solution()
{
  ofstream some_file;
  char filename[100];

  // Output the time
  cout << "Time is now " << time_pt()->time() << std::endl;

  // Number of plot points
  unsigned npts;
  npts = 5;

  sprintf(
    filename, "%s/soln%i.dat", Doc_info.directory().c_str(), Doc_info.number());
  some_file.open(filename);
  this->My_mesh_pt->output(some_file, npts);
  some_file.close();

  ofstream timefile;
  timefile.open(GlobalSimSettings::result_folder + "/time.dat",
                std::ofstream::out | std::ofstream::app);
  timefile << time_pt()->time() << "\n";
  timefile.close();

  // Need a measure of how close the fronts are at late times.
  unsigned nnod = My_mesh_pt->nnode();
  // First find the inner front
  Vector<unsigned> front_index(2, 0);
  unsigned k = 0;
  // for (unsigned n = 11; n < nnod - 11; n++)
  for (unsigned n = nnod - 12; n > 10; --n)
  {
    // Get the concentration at node: average over nearest 20 neighbours
    double phi_local = 0.0;
    for (unsigned n2 = n - 10; n2 < n + 10; n2++)
    {
      phi_local += My_mesh_pt->node_pt(n2)->value(2);
    }
    phi_local = phi_local / 20.0;

    double phi_thres = 0.95 * GlobalPhysicalVariables::phi_c;
    if (phi_local > phi_thres)
    {
      // Look at adjacent nodes: do the same averaging
      double phi_plus = 0.0;
      double phi_minus = 0.0;
      for (unsigned n2 = n - 10; n2 < n + 10; n2++)
      {
        phi_plus += My_mesh_pt->node_pt(n2 + 1)->value(2);
        phi_minus += My_mesh_pt->node_pt(n2 - 1)->value(2);
      }
      phi_plus = phi_plus / 20.0;
      phi_minus = phi_minus / 20.0;

      // If either of these are below the threshold, then we have a front.
      // Record it.
      if ((phi_plus < phi_thres || phi_minus < phi_thres) && (k < 2))
      {
        front_index[k] = n;
        ++k;
      }
    }
  }

  // Doc the fronts
  ofstream fronts_file;
  char fronts_filename[100];
  sprintf(fronts_filename,
          "%s/fronts%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  fronts_file.open(fronts_filename);
  fronts_file << My_mesh_pt->node_pt(front_index[0])->x(0) << std::endl;
  fronts_file << My_mesh_pt->node_pt(front_index[1])->x(0) << std::endl;
  fronts_file.close();


  // If there are two distinct fronts and they are close enough, reduce the
  // Peclet number at late times
  if ((front_index[0] != front_index[1]) && (front_index[0] != 0) &&
      (front_index[1] != 0))
  {
    double delta = std::fabs(My_mesh_pt->node_pt(front_index[0])->x(0) -
                             My_mesh_pt->node_pt(front_index[1])->x(0));
    if (delta < 0.5)
    {
      GlobalPhysicalVariables::Pe = 50.0;
    }
  }

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc


//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template<class ELEMENT>
double ThinFilmProblem<ELEMENT>::global_temporal_error_norm()
{
  double global_error = 0.0;

  // Find out how many nodes there are in the problem
  unsigned n_node = My_mesh_pt->nnode();

  // Loop over the nodes and calculate the errors in the values
  for (unsigned i = 0; i < n_node; i++)
  {
    // Get error in solution: Difference between predicted and actual
    // value for nodal value 0
    double error =
      My_mesh_pt->node_pt(i)->time_stepper_pt()->temporal_error_in_value(
        My_mesh_pt->node_pt(i), 0);

    // Add the square of the individual error to the global error
    global_error += error * error;
  }

  // Now the global error must be divided by the number of nodes
  global_error /= double(n_node);

  // Return the square root of the error
  return sqrt(global_error);

} // end of global_temporal_error_norm


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char** argv)
{
  //--------Store command line arguments-----------------------//
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--peclet", &GlobalPhysicalVariables::Pe, "Peclet number");

  CommandLineArgs::specify_command_line_flag(
    "--nu", &GlobalPhysicalVariables::nu, "Inverse pore size");

  CommandLineArgs::specify_command_line_flag(
    "--capillary", &GlobalPhysicalVariables::Ca, "Capillary number");

  CommandLineArgs::specify_command_line_flag("--phi_initial",
                                             &GlobalPhysicalVariables::phi_0,
                                             "Initial volume fraction");

  CommandLineArgs::specify_command_line_flag(
    "--folder", &GlobalSimSettings::result_folder, "Result folder");

  CommandLineArgs::specify_command_line_flag(
    "--t_max", &GlobalSimSettings::t_max, "Max time");

  CommandLineArgs::specify_command_line_flag(
    "--evaporation_flag",
    &GlobalPhysicalVariables::is_diffusive,
    "Evaporation flag");

  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  //-----------------------------------------------------------//
  // Write all parameters to a file
  //-----------------------------------------------------------//
  ofstream parameter_file(GlobalSimSettings::result_folder + "/parameters.dat");
  parameter_file << "Pe: " << GlobalPhysicalVariables::Pe << std::endl;
  parameter_file << "nu: " << GlobalPhysicalVariables::nu << std::endl;
  parameter_file << "Ca: " << GlobalPhysicalVariables::Ca << std::endl;
  parameter_file << "phi_0: " << GlobalPhysicalVariables::phi_0 << std::endl;
  parameter_file << "Max permitted error: "
                 << GlobalSimSettings::max_permitted_error << std::endl;
  parameter_file << "Min permitted error: "
                 << GlobalSimSettings::min_permitted_error << std::endl;
  if (GlobalPhysicalVariables::is_diffusive)
  {
    parameter_file << "Evaporation mode: diffusive" << std::endl;
  }
  else
  {
    parameter_file << "Evaporation mode: kinetic" << std::endl;
  }
  parameter_file << std::endl;
  parameter_file << "Transition function parameters:" << std::endl;
  parameter_file << "delta: " << GlobalPhysicalVariables::delta << std::endl;
  parameter_file << "phi_tilde: " << GlobalPhysicalVariables::phi_tilde
                 << std::endl;
  parameter_file.close();
  //-----------------------------------------------------------//

  // Clear timefile before start
  ofstream timefile(GlobalSimSettings::result_folder + "/time.dat");

  // Create problem
  ThinFilmProblem<RefineableQThinFilmBrinkmanElement<1, 3>> problem;

  // Set the initial condition *after* uniform refinement
  problem.set_initial_condition();

  // Doc initial solution
  problem.doc_solution();

  // Set the initial timestep
  double dt = GlobalSimSettings::dt_comp;

  // Initial target time
  double t_target = 0.05;

  bool increment_boolean = false;
  double dt_max = 0.0095;

  // Timestepping loop:
  unsigned k = 0;
  while (problem.time_pt()->time() < GlobalSimSettings::t_max)
  {
    double dt_new =
      problem.adaptive_unsteady_newton_solve(dt, GlobalSimSettings::epsilon_t);

    dt_new = std::min(dt_new, dt_max);

    // The actual timestep taken:
    double dt_taken = problem.time_pt()->time() - problem.time_pt()->time(1);

    // If we are due to increment and the actual timestep was dt, increment
    // t_target by 0.1. Also set the full output boolean to true.
    if (increment_boolean && fabs(dt_taken - dt) < 1.0e-6)
    {
      t_target += 0.05;
      increment_boolean = false;
    }

    // If the next trial timestep exceeds the next desired output time, shorten
    // accordingly
    if (problem.time_pt()->time() + dt_new > t_target)
    {
      dt = t_target - problem.time_pt()->time();
      increment_boolean = true;
    }
    else
    {
      dt = dt_new;
    }

    if (k > 9)
    {
      problem.adapt();
    }
    ++k;
    
    // Doc the solution
    problem.doc_solution();
  }

} // End of main
