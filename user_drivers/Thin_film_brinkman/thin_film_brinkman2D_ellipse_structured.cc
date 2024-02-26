//====================================================================
//  Driver for a 2D drying drop problem on a structured mesh
//====================================================================

#include <fenv.h>

// Generic routines
#include "generic.h"

// The equations
#include "thin_film_brinkman.h"

// Elliptic integrals
#include "elliptic_integral.h"

// The meshes
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;
using namespace oomph;

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
  double phi_0 = phi_c / 25.0;

  /// Peclet number (specified as a command line flag)
  double Pe;

  /// Scaled inverse pore size
  double nu = 1.0e3;

  /// Axisymmetric geometry?
  bool axisym_flag = false;

  /// Flag for kinetic/diffusive evaporation.
  unsigned is_diffusive = 0;

  // Parameters for the transition function
  double phi_tilde = 0.7 * phi_c;
  double delta = 0.001;
  double alpha;
  double beta;

  // The "dry-out" time in the absence of jamming
  double t_f;

  /// Ellipse major and minor axes
  double a = 1.0;
  double b = 1.0;

  /// Pinning height at contact line
  double pinning_height = 0.0;

  /// Initial height of the drop
  double h0;

  /// Evaporation flux function
  void get_evaporation_flux_fct(const double& time,
                                const Vector<double>& x,
                                double& evap)
  {
    if (is_diffusive)
    {
      // First get the eccentricity and pi
      double pi = MathematicalConstants::Pi;
      double ecc = pow(1.0 - b * b / (a * a), 0.5);

      evap = (pi / (4.0 * b * elliptic_fk(ecc))) *
             pow(1.0 - pow(x[0] / a, 2.0) - pow(x[1] / b, 2.0), -0.5);
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

    Q = -1.0 * (mu_KD_inv * pow(h, 3.0) / (3.0 * Ca) + h / (Ca * pow(nu, 2.0)));
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

    dQdh = -1.0 * (pow(h, 2.0) * mu_KD_inv / Ca + 1.0 / (pow(nu, 2.0) * Ca));
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

    dQdphi = -1.0 * pow(h, 3.0) * mu_KD_inv_deriv / (3 * Ca);
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

    Q = -1.0 * transition_fct * mu_KD_inv * pow(h, 3.0) / (3.0 * Ca);
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


    dQdh = -1.0 * transition_fct * mu_KD_inv * pow(h, 2.0) / Ca;
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
      -1.0 *
      (transition_fct_deriv * mu_KD_inv + transition_fct * mu_KD_inv_deriv) *
      pow(h, 3.0) / (3.0 * Ca);
  }

} // namespace GlobalPhysicalVariables


//======start_of_GlobalSimSettings====================================
/// Namespace for global simulation settings
//====================================================================
namespace GlobalSimSettings
{
  /// Result Folder
  string result_folder = "RESLT_2D_ellipse";

  /// Are we picking up from a restart?
  unsigned restart = 0;

  // Min. refinement level: default is zero
  unsigned min_refinement_level = 3;

  // Max. refinement level: default is five
  unsigned max_refinement_level = 11;

  // Maximum permitted error target (elements refined if error is above this)
  double max_permitted_error = 5.0e-4;

  // Minimum permitted error target (elements unrefined if error is below this)
  double min_permitted_error = 1.0e-4;

  // Target error for adaptive timestepping
  double epsilon_t = 1.0e-3;

  // How frequently the mesh is rebuilt
  unsigned adapt_number = 3;

  // Number of output timesteps
  unsigned Nsteps = 100;

  // Maximum simulation time
  double t_max = 0.6;

  // Computational timestep
  double dt_comp = 1.0e-5;

} // namespace GlobalSimSettings


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class StructuredThinFilmProblem : public virtual Problem
{
public:
  /// Constructor
  StructuredThinFilmProblem();

  /// Destructor
  ~StructuredThinFilmProblem()
  {
    delete outer_boundary_ellipse_pt;
    delete My_mesh_pt;
    delete error_estimator_pt;
    delete this->time_stepper_pt();
  };

  /// Actions before adapt (empty)
  void actions_before_adapt() {}

  /// Update the problem specs after timestep (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep (empty)
  void actions_before_implicit_timestep() {}

  /// Set initial condition
  void set_initial_condition();

  /// Actions after adapt:
  /// Setup the problem again
  void actions_after_adapt()
  {
    complete_problem_setup();
  }

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}

  /// Doc the solution: boolean flags for full output and dumping the
  /// solution in case of a crash
  void doc_solution(bool full_output);

  /// Global error norm for adaptive time-stepping
  double global_temporal_error_norm();

  /// Dump problem to disk to allow for restart.
  void dump_it(ofstream& dump_file);

  /// Read problem for restart from specified restart file.
  void restart(ifstream& restart_file);

  /// Wrapper for adaptive timestepping
  void custom_adaptive_unsteady_newton_solve(const double& dt_desired,
                                             const double& epsilon,
                                             double& dt_new,
                                             bool& success);

  /// Access function for the mesh:
  RefineableQuarterCircleSectorMesh<ELEMENT>* my_mesh_pt()
  {
    return My_mesh_pt;
  }

  void smooth_spurious_nodes();

private:
  /// Doc info object for labeling output
  DocInfo Doc_info;

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of  all elements
  void complete_problem_setup();

  /// Doc uniform mesh during full output
  void doc_uniform_mesh();

  /// Pointers to specific mesh
  RefineableQuarterCircleSectorMesh<ELEMENT>* My_mesh_pt;

  /// A uniformly refined mesh purely for output
  RefineableQuarterCircleSectorMesh<ELEMENT>* Output_mesh_pt;

  /// Pointers for mesh construction
  Ellipse* outer_boundary_ellipse_pt;

  // Pointer to error estimator
  Z2ErrorEstimator* error_estimator_pt;
}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor
//========================================================================
template<class ELEMENT>
StructuredThinFilmProblem<ELEMENT>::StructuredThinFilmProblem()
{
  // Output directory
  Doc_info.set_directory(GlobalSimSettings::result_folder);

  // Output number
  Doc_info.number() = 0;

  // Allocate the timestepper (this constructs the time object as well)
  // Boolean flag enables temporal adaptivity
  add_time_stepper_pt(new BDF<2>(true));

#ifdef OOMPH_HAS_MUMPS
  // Use mumps if available
  linear_solver_pt() = new MumpsSolver;
#endif

  // Don't reject timesteps above tolerance
  Keep_temporal_error_below_tolerance = false;

  // Intrinsic coordinate along GeomObject
  Vector<double> zeta(1);

  // Position vector on GeomObject
  Vector<double> posn(2);

  double xi_lo = 0.0;
  double xi_hi = MathematicalConstants::Pi / 2.0;
  double fract_mid = 0.5;

  // Ellipse defining the closing cap
  Ellipse* outer_boundary_ellipse_pt =
    new Ellipse(GlobalPhysicalVariables::a, GlobalPhysicalVariables::b);

  // Create the mesh
  My_mesh_pt = new RefineableQuarterCircleSectorMesh<ELEMENT>(
    outer_boundary_ellipse_pt, xi_lo, fract_mid, xi_hi, time_stepper_pt());

  // Also create the output mesh
  Output_mesh_pt = new RefineableQuarterCircleSectorMesh<ELEMENT>(
    outer_boundary_ellipse_pt, xi_lo, fract_mid, xi_hi, time_stepper_pt());

  // Refine output mesh
  for (unsigned n = 0; n < 5; n++)
  {
    Output_mesh_pt->refine_uniformly();
  }

  // Bulk mesh is the problem's only mesh
  Problem::mesh_pt() = My_mesh_pt;

  // Set error estimator for bulk mesh
  error_estimator_pt = new Z2ErrorEstimator;
  My_mesh_pt->spatial_error_estimator_pt() = error_estimator_pt;

  // Set the error targets
  My_mesh_pt->max_permitted_error() = GlobalSimSettings::max_permitted_error;
  My_mesh_pt->min_permitted_error() = GlobalSimSettings::min_permitted_error;

  // Max. refinement level
  My_mesh_pt->min_refinement_level() = GlobalSimSettings::min_refinement_level;
  My_mesh_pt->max_refinement_level() = GlobalSimSettings::max_refinement_level;

  // Set boundary condition and complete the build of all elements
  complete_problem_setup();

  // Refine uniformly twice to begin with
  if (!GlobalSimSettings::restart)
  {
    for (unsigned n = 0; n < GlobalSimSettings::min_refinement_level; n++)
    {
      refine_uniformly();
    }
  }

  // Setup equation numbering scheme
  oomph_info << "Number of equations: " << this->assign_eqn_numbers()
             << std::endl;

} // end_of_constructor


//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::complete_problem_setup()
{
  // Set the boundary conditions for problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  unsigned num_nod = My_mesh_pt->nboundary_node(1);
  for (unsigned inod = 0; inod < num_nod; inod++)
  {
    // Get node
    Node* nod_pt = My_mesh_pt->boundary_node_pt(1, inod);

    // Pin the free surface height only. Leaving the other variable
    // un-pinned imposes a Neumann condition on h
    nod_pt->pin(0);
  }

  // Complete the build of all elements so they are fully functional
  unsigned n_element = My_mesh_pt->nelement();

  for (unsigned e = 0; e < n_element; e++)
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

    // Set the solute mobility derivative dQdh pointer
    el_pt->solute_mobility_fct_deriv_h_pt() =
      &GlobalPhysicalVariables::get_partial_solute_mobility_partial_h;

    // Set the solute mobility derivative dQdphi pointer
    el_pt->solute_mobility_fct_deriv_c_pt() =
      &GlobalPhysicalVariables::get_partial_solute_mobility_partial_phi;

    // Set the Peclet number pointer
    el_pt->peclet_pt() = &GlobalPhysicalVariables::Pe;

    // Set the dry-out time
    el_pt->t_f_pt() = &GlobalPhysicalVariables::t_f;

    // Set the axisymmetry flag
    el_pt->axisymmetry_flag_pt() = &GlobalPhysicalVariables::axisym_flag;
  }

  // Re-apply Dirichlet boundary conditions
  apply_boundary_conditions();
}


//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::apply_boundary_conditions()
{
  // Loop over all boundary nodes
  unsigned num_nod = this->My_mesh_pt->nboundary_node(1);
  for (unsigned inod = 0; inod < num_nod; inod++)
  {
    // Get node
    Node* nod_pt = this->My_mesh_pt->boundary_node_pt(1, inod);

    // Fix the free surface height to the pinning height
    nod_pt->set_value(0, GlobalPhysicalVariables::pinning_height);
  }
} // end set bc


//==start_of_set_initial_condition========================================
/// Set initial conditions: Set all nodal velocities to zero and
/// initialise the previous velocities and nodal positions to correspond
/// to an impulsive start
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::set_initial_condition()
{
  // Pointer to restart file
  ifstream* restart_file_pt = 0;

  // Restart?
  //---------
  // Restart file specified via command line [all programs have at least
  // a single command line argument: their name. Ignore this here.]
  if (GlobalSimSettings::restart)
  {
    // Open restart file
    restart_file_pt = new ifstream(
      GlobalSimSettings::result_folder + "/restart.dat", ios_base::in);
    if (restart_file_pt != 0)
    {
      oomph_info << "Have set initial condition from a restart." << std::endl;
    }
    else
    {
      std::ostringstream error_stream;
      error_stream << "Restart file not found." << std::endl;

      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }
  else
  {
    oomph_info << "No restart -- setting IC as normal" << std::endl;
  }

  // Read restart data:
  //-------------------
  if (restart_file_pt != 0)
  {
    // Read the problem data from the restart file
    restart(*restart_file_pt);
  }
  // Assign initial condition as normal
  //---------------------------------------------
  else
  {
    // Initialise the timestep
    initialise_dt(GlobalSimSettings::dt_comp);

    // Find number of nodes in mesh
    unsigned num_nod = My_mesh_pt->nnode();

    Vector<double> x(2);

    Node* nod_pt;

    // Loop over the nodes to set initial guess everywhere
    for (unsigned jnod = 0; jnod < num_nod; jnod++)
    {
      // Get node
      nod_pt = My_mesh_pt->node_pt(jnod);

      x[0] = nod_pt->x(0);
      x[1] = nod_pt->x(1);

      // Initial solute volume fraction is phi_0
      mesh_pt()->node_pt(jnod)->set_value(2, GlobalPhysicalVariables::phi_0);

      // Initial height is the poisson solution, with apex height h0
      mesh_pt()->node_pt(jnod)->set_value(
        0,
        GlobalPhysicalVariables::h0 *
          (1.0 - pow(x[0] / GlobalPhysicalVariables::a, 2.0) -
           pow(x[1] / GlobalPhysicalVariables::b, 2.0)));
    }

    // History values correspond to an impulsive start
    assign_initial_values_impulsive();
  }
} // End of set_initial_condition


//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
template<class ELEMENT>
double StructuredThinFilmProblem<ELEMENT>::global_temporal_error_norm()
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


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::doc_solution(bool full_output)
{
  // Always output data on the long and short axes
  ofstream some_fileX;
  ofstream some_fileY;
  char filenameX[100];
  char filenameY[100];

  // Output the time
  oomph_info << "Time is now " << time_pt()->time() << std::endl;

  sprintf(filenameX,
          "%s/long_axis%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  sprintf(filenameY,
          "%s/short_axis%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_fileX.open(filenameX);
  some_fileY.open(filenameY);

  // Create temporary vectors which hold the coordinates and concentrations
  unsigned n_nod = My_mesh_pt->nboundary_node(0);
  vector<int> x_indices(n_nod);
  vector<double> x_coords(n_nod);
  // Assign all values
  for (unsigned n = 0; n < n_nod; n++)
  {
    x_indices[n] = n;
    x_coords[n] = My_mesh_pt->boundary_node_pt(0, n)->x(0);
  }
  // Now sort
  sort(x_indices.begin(),
       x_indices.end(),
       [&](int i, int j) { return x_coords[i] < x_coords[j]; });

  // Do the same on boundary 2
  n_nod = My_mesh_pt->nboundary_node(2);
  vector<int> y_indices(n_nod);
  vector<double> y_coords(n_nod);
  for (unsigned n = 0; n < n_nod; n++)
  {
    y_indices[n] = n;
    y_coords[n] = My_mesh_pt->boundary_node_pt(2, n)->x(1);
  }
  sort(y_indices.begin(),
       y_indices.end(),
       [&](int i, int j) { return y_coords[i] < y_coords[j]; });

  // Output everything on the long and short axes
  for (unsigned b = 0; b < 3; b++)
  {
    n_nod = My_mesh_pt->nboundary_node(b);
    for (unsigned n = 0; n < n_nod; n++)
    {
      if (b == 0)
      {
        unsigned tmp_idx = x_indices[n];
        some_fileX << My_mesh_pt->boundary_node_pt(b, tmp_idx)->x(0) << " ";
        for (unsigned a = 0; a < 3; a++)
        {
          some_fileX << My_mesh_pt->boundary_node_pt(b, tmp_idx)->value(a)
                     << " ";
        }
        some_fileX << std::endl;
      }
      if (b == 2)
      {
        unsigned tmp_idx = y_indices[n];
        some_fileY << My_mesh_pt->boundary_node_pt(b, tmp_idx)->x(1) << " ";
        for (unsigned a = 0; a < 3; a++)
        {
          some_fileY << My_mesh_pt->boundary_node_pt(b, tmp_idx)->value(a)
                     << " ";
        }
        some_fileY << std::endl;
      }
    }
  }
  some_fileX.close();
  some_fileY.close();

  // Need a measure of how close the fronts are at late times.
  n_nod = My_mesh_pt->nboundary_node(0);
  // First find the inner front
  Vector<unsigned> front_index(2, 0.0);
  unsigned k = 0;
  for (unsigned n = 4; n < n_nod - 4; n++)
  {
    unsigned tmp_idx = x_indices[n];

    // Get the concentration at node: average over nearest 6 neighbours
    double phi_local = 0.0;
    for (unsigned n2 = n - 3; n2 < n + 3; n2++)
    {
      unsigned tmp_idx2 = x_indices[n2];
      phi_local += My_mesh_pt->boundary_node_pt(0, tmp_idx2)->value(2);
    }
    phi_local = phi_local / 6.0;

    double phi_thres = 0.95 * GlobalPhysicalVariables::phi_c;
    if (phi_local > phi_thres)
    {
      // Look at adjacent nodes: do the same averaging
      double phi_plus = 0.0;
      double phi_minus = 0.0;
      for (unsigned n2 = n - 3; n2 < n + 3; n2++)
      {
        unsigned tmp_idx_p = x_indices[n2 + 1];
        unsigned tmp_idx_m = x_indices[n2 - 1];
        phi_plus += My_mesh_pt->boundary_node_pt(0, tmp_idx_p)->value(2);
        phi_minus += My_mesh_pt->boundary_node_pt(0, tmp_idx_m)->value(2);
      }
      phi_plus = phi_plus / 6.0;
      phi_minus = phi_minus / 6.0;

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
  fronts_file
    << My_mesh_pt->boundary_node_pt(0, x_indices[front_index[0]])->x(0)
    << std::endl;
  fronts_file
    << My_mesh_pt->boundary_node_pt(0, x_indices[front_index[1]])->x(0)
    << std::endl;
  fronts_file.close();

  // If there are two distinct fronts and they are close enough, reduce the
  // Peclet number at late times
  if ((front_index[0] != front_index[1]) && (front_index[0] != 0) &&
      (front_index[1] != 0))
  {
    double delta = std::fabs(
      My_mesh_pt->boundary_node_pt(0, x_indices[front_index[0]])->x(0) -
      My_mesh_pt->boundary_node_pt(0, x_indices[front_index[1]])->x(0));
    if (delta < 0.6 * GlobalPhysicalVariables::a)
    {
      GlobalPhysicalVariables::Pe = 20.0;
    }
  }

  // Output current time
  ofstream timefile;
  timefile.open(GlobalSimSettings::result_folder + "/time.dat",
                std::ofstream::out | std::ofstream::app);
  timefile << time_pt()->time() << "\n";
  timefile.close();

  // If we're due, do a full output
  if (full_output)
  {
    ofstream some_file;
    char filename[100];

    // Output the time
    oomph_info << "Doing a full output and (over)writing restart file."
               << std::endl;

    // Number of plot points
    unsigned npts;
    npts = 3;

    sprintf(filename,
            "%s/full_soln%i.dat",
            Doc_info.directory().c_str(),
            Doc_info.number());

    some_file.open(filename);
    this->My_mesh_pt->output(some_file, npts);
    some_file.close();

    // Also write to restart file in case of crash. NB this will overwrite
    // the restart file every time.
    sprintf(filename, "%s/restart.dat", Doc_info.directory().c_str());
    some_file.open(filename);
    dump_it(some_file);
    some_file.close();

    // Doc the uniform mesh for processing output
    doc_uniform_mesh();
  }

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc


//=====start_of_dump_it===================================================
/// Dump the solution to disk to allow for restart
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::dump_it(ofstream& dump_file)
{
  // Call generic dump()
  Problem::dump(dump_file);
} // end of dump_it


//=================start_of_restart=======================================
/// Read solution from disk for restart
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::restart(ifstream& restart_file)
{
  // Read the generic problem data from restart file
  Problem::read(restart_file);
} // end of restart


//===============start_of_custom_adaptive_timestepper=====================
/// A slight variation of adaptive_unsteady_newton_solve. Does a doc when
/// the solver fails
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::custom_adaptive_unsteady_newton_solve(
  const double& dt_desired,
  const double& epsilon,
  double& dt_new,
  bool& success)
{
  // Assume we will succeed
  success = true;

  // First, we need to backup the existing dofs, in case the timestep is
  // rejected

  // Find total number of dofs on current processor
  unsigned n_dof_local = dof_distribution_pt()->nrow_local();

  // Now set up a Vector to hold current values
  Vector<double> dofs_current(n_dof_local);

  // Load values into dofs_current
  for (unsigned i = 0; i < n_dof_local; i++) dofs_current[i] = dof(i);

  // Store the time
  double time_current = time_pt()->time();

  // Flag to detect whether the timestep has been rejected or not
  bool reject_timestep = 0;

  // Flag to detect whether any of the timesteppers are adaptive
  unsigned adaptive_flag = 0;

  // The value of the actual timestep, by default the same as desired timestep
  double dt_actual = dt_desired;

  // Find out whether any of the timesteppers are adaptive
  unsigned n_time_steppers = ntime_stepper();
  for (unsigned i = 0; i < n_time_steppers; i++)
  {
    if (time_stepper_pt(i)->adaptive_flag())
    {
      adaptive_flag = 1;
      break;
    }
  }

  // Shift the time_values
  shift_time_values();

  // This loop surrounds the adaptive time-stepping and will not be broken
  // until a timestep is accepted
  do
  {
    // Initially we assume that this step will succeed and that this dt
    // value is ok.
    reject_timestep = 0;
    double dt_rescaling_factor = 1.0;

    // Set the new time and value of dt
    time_pt()->time() += dt_actual;
    time_pt()->dt() = dt_actual;

    // Loop over all timesteppers and set the weights and predictor weights
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      // If the time_stepper is non-adaptive, this will be zero
      time_stepper_pt(i)->set_predictor_weights();
      time_stepper_pt(i)->set_weights();
    }

    // Now calculate the predicted values for the all data and all positions
    calculate_predictions();

    // Run the individual timesteppers actions before timestep. These need to
    // be before the problem's actions_before_implicit_timestep so that the
    // boundary conditions are set consistently.
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->actions_before_timestep(this);
    }

    // Do any updates/boundary conditions changes here
    actions_before_implicit_timestep();

    // Attempt to solve the non-linear system
    try
    {
      // Solve the non-linear problem at this timestep
      newton_solve();
    }
    // Catch any exceptions thrown
    catch (NewtonSolverError& error)
    {
      // If it's a solver error then die
      if (error.linear_solver_error || Time_adaptive_newton_crash_on_solve_fail)
      {
        std::string error_message = "USER-DEFINED ERROR IN NEWTON SOLVER\n";
        error_message += "ERROR IN THE LINEAR SOLVER\n";

        // Die
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Reject the timestep, if we have an exception
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;

        // Half the time step
        dt_rescaling_factor = Timestep_reduction_factor_after_nonconvergence;
      }
    }

    // Run the individual timesteppers actions, these need to be before the
    // problem's actions_after_implicit_timestep so that the time step is
    // finished before the problem does any auxiliary calculations (e.g. in
    // semi-implicit micromagnetics the calculation of magnetostatic field).
    for (unsigned i = 0; i < n_time_steppers; i++)
    {
      time_stepper_pt(i)->actions_after_timestep(this);
    }

    // Update anything that needs updating after the timestep
    actions_after_implicit_timestep();

    // If we have an adapative timestepper (and we haven't already failed)
    // then calculate the error estimate and rescaling factor.
    if (adaptive_flag && !reject_timestep)
    {
      // Once timestep has been accepted can do fancy error processing
      // Set the error weights
      for (unsigned i = 0; i < n_time_steppers; i++)
      {
        time_stepper_pt(i)->set_error_weights();
      }

      // Get a global error norm to use in adaptivity (as specified by the
      // problem sub-class writer). Prevent a divide by zero if the solution
      // gives very close to zero error. Error norm should never be negative
      // but use absolute value just in case.
      double error = std::max(std::abs(global_temporal_error_norm()), 1e-12);

      // Calculate the scaling  factor
      dt_rescaling_factor =
        std::pow((epsilon / error), (1.0 / (1.0 + time_stepper_pt()->order())));

      oomph_info << "Timestep scaling factor is  " << dt_rescaling_factor
                 << std::endl;
      oomph_info << "Estimated timestepping error is " << error << std::endl;


      // Do we have to do it again?
      if (error > epsilon)
      {
        oomph_info << "Estimated timestepping error " << error
                   << " exceeds tolerance " << epsilon << " ";
        if (Keep_temporal_error_below_tolerance)
        {
          oomph_info << " --> rejecting timestep." << std::endl;
          reject_timestep = 1;
        }
        else
        {
          oomph_info << " ...but we're not rejecting the timestep" << std::endl;
        }
        oomph_info
          << "Note: This behaviour can be adjusted by changing the protected "
          << "boolean" << std::endl
          << std::endl
          << "    Problem::Keep_temporal_error_below_tolerance" << std::endl;
      }


    } // End of if adaptive flag


    // Calculate the next time step size and check it's ok
    // ============================================================

    // Calculate the possible next time step, if no error conditions
    // trigger.
    double new_dt_candidate = dt_rescaling_factor * dt_actual;
    double custom_time_scaling_factor = 2.0;

    // Check that the scaling factor is within the allowed range
    if (dt_rescaling_factor > custom_time_scaling_factor)
    {
      oomph_info << "Tried to increase dt by the ratio " << dt_rescaling_factor
                 << " which is above the maximum ("
                 << custom_time_scaling_factor
                 << "). Attempting to increase by the maximum ratio instead."
                 << std::endl;
      new_dt_candidate = custom_time_scaling_factor * dt_actual;
    }
    // If we have already rejected the timestep then don't do this check
    // because DTSF will definitely be too small.
    else if ((!reject_timestep) && (dt_rescaling_factor <= DTSF_min_decrease))
    {
      // Handle this special case where we want to continue anyway (usually
      // Minimum_dt_but_still_proceed = -1 so this has no effect).
      if (new_dt_candidate < Minimum_dt_but_still_proceed)
      {
        oomph_info << "Warning: Adaptation of timestep to ensure satisfaction\n"
                   << "         of error bounds during adaptive timestepping\n"
                   << "         would lower dt below \n"
                   << "         Problem::Minimum_dt_but_still_proceed="
                   << Minimum_dt_but_still_proceed << "\n"
                   << "         ---> We're continuing with present timestep.\n"
                   << std::endl;
        dt_rescaling_factor = 1.0;
        // ??ds shouldn't we set new_dt_candidate =
        // Minimum_dt_but_still_proceed here, rather than not changing dt at
        // all?
      }
      else
      {
        // Otherwise reject
        oomph_info << "Timestep would decrease by " << dt_rescaling_factor
                   << " which is less than the minimum scaling factor "
                   << DTSF_min_decrease << std::endl;
        oomph_info << "TIMESTEP REJECTED" << std::endl;
        reject_timestep = 1;
      }
    }

    // Now check that the new dt is within the allowed range
    if (new_dt_candidate > Maximum_dt)
    {
      oomph_info << "Tried to increase dt to " << new_dt_candidate
                 << " which is above the maximum (" << Maximum_dt
                 << "). I increased it to the maximum value instead.";
      dt_actual = Maximum_dt;
    }
    else if (new_dt_candidate < Minimum_dt)
    {
      oomph_info << "MINIMUM TIMESTEP HAS BEEN REACHED: BREAKING OUT AND DOING "
                    "A FINAL FULL DOC."
                 << std::endl;
      doc_solution(true);
      success = false;
      break;
    }
    else
    {
      dt_actual = new_dt_candidate;
    }


    actions_after_implicit_timestep_and_error_estimation();


    // If we are rejecting this attempt then revert the dofs etc.
    if (reject_timestep)
    {
      // Reset the time
      time_pt()->time() = time_current;

      // Reload the dofs
      unsigned ni = dofs_current.size();
      for (unsigned i = 0; i < ni; i++)
      {
        dof(i) = dofs_current[i];
      }

#ifdef OOMPH_HAS_MPI
      // Synchronise the solution on different processors (on each submesh)
      this->synchronise_all_dofs();
#endif

      // Call all "after" actions, e.g. to handle mesh updates
      actions_after_newton_step();
      actions_before_newton_convergence_check();
      actions_after_newton_solve();
      actions_after_implicit_timestep();
      actions_after_implicit_timestep_and_error_estimation();
    }

  }
  // Keep this loop going until we accept the timestep
  while (reject_timestep);

  // Once the timestep has been accepted, return the time step that should be
  // used next time.
  dt_new = dt_actual;
}

//==start_of_doc_uniform_mesh=============================================
/// Project onto a uniform mesh via the locate_zeta function and output
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::doc_uniform_mesh()
{
#ifdef OOMPH_HAS_CGAL

  // Make cgal-based bin
  CGALSamplePointContainerParameters cgal_params(My_mesh_pt);
  MeshAsGeomObject* mesh_geom_obj_pt = new MeshAsGeomObject(&cgal_params);

  // Open output file
  ofstream some_file;
  char filename[100];
  sprintf(filename,
          "%s/full_soln_uniform%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());

  some_file.open(filename);

  // Number of 1D output points per element
  unsigned n_pts = 5;

  // Containers for local and global coordinates of output mesh
  Vector<double> s(2, 0.0);
  Vector<double> x(2, 0.0);

  // Container for local coordinate in computational mesh
  Vector<double> s_comp(2, 0.0);

  // Loop over elements in the output mesh
  unsigned n_el = Output_mesh_pt->nelement();
  for (unsigned n = 0; n < n_el; n++)
  {
    // Cast to the element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Output_mesh_pt->element_pt(n));

    // Loop twice up to n_pts, one loop for each dimension
    for (unsigned i = 0; i < n_pts; i++)
    {
      for (unsigned j = 0; j < n_pts; j++)
      {
        // The local coordinate s
        s[0] = -1.0 + 2.0 * i / (n_pts - 1);
        s[1] = -1.0 + 2.0 * j / (n_pts - 1);

        // The global coords
        el_pt->interpolated_x(s, x);

        // Allocate temporary geomobject
        GeomObject* geom_object_pt = 0;

        // Locate the coordinate in computational mesh
        dynamic_cast<CGALSamplePointContainer*>(
          mesh_geom_obj_pt->sample_point_container_pt())
          ->locate_zeta(x, geom_object_pt, s_comp);

        if (geom_object_pt)
        {
          // Up-cast to element
          ELEMENT* el2_pt = dynamic_cast<ELEMENT*>(geom_object_pt);

          // Output all fields
          some_file << x[0] << " " << x[1] << " ";
          for (unsigned ifield = 0; ifield < 3; ifield++)
          {
            some_file << el2_pt->interpolated_u_thin_film_brinkman(s_comp,
                                                                   ifield)
                      << " ";
          }
          
          some_file << std::endl;
        }
      }
    }
  }

  // Delete temporary mesh geom object
  delete mesh_geom_obj_pt;

  // Close output file
  some_file.close();

#endif
}


//==start_of_smooth_spurious nodes========================================
/// Nodes with spurious phi values are set to a local nodal average
//========================================================================
template<class ELEMENT>
void StructuredThinFilmProblem<ELEMENT>::smooth_spurious_nodes()
{
  // Any nodes with an instability: set them to nearest neighbour value
  unsigned n_nod = My_mesh_pt->nnode();
  for (unsigned n = 0; n < n_nod; n++)
  {
    Node* nod_pt = My_mesh_pt->node_pt(n);
    if (nod_pt->value(2) > 0.642)
    {
      Vector<unsigned> nearest_nod_idxs(3, 1000000000);
      for (unsigned k = 0; k < 3; k++)
      {
        double r = 1.0e5;

        for (unsigned n2 = 0; n2 < n_nod; n2++)
        {
          if ((n2 != n) && (n2 != nearest_nod_idxs[0]) &&
              (n2 != nearest_nod_idxs[1]) && (n2 != nearest_nod_idxs[2]))
          {
            Node* nod2_pt = My_mesh_pt->node_pt(n2);
            double r_tmp = pow(pow(nod_pt->x(0) - nod2_pt->x(0), 2.0) +
                                 pow(nod_pt->x(1) - nod2_pt->x(1), 2.0),
                               0.5);
            if (r_tmp < r)
            {
              r = r_tmp;
              nearest_nod_idxs[k] = n2;
            }
          }
        }
      }

      double phi_val = 0.0;
      for (unsigned k = 0; k < 3; k++)
      {
        phi_val += My_mesh_pt->node_pt(nearest_nod_idxs[k])->value(2);
      }
      phi_val = phi_val / 3.0;

      nod_pt->set_value(2, phi_val);
    }
  }
}


//=======start_of_main========================================
//============================================================
int main(int argc, char** argv)
{
  // Start time of the run
  time_t start_time = time(NULL);

#ifdef OOMPH_HAS_MPI

  // Initialise MPI
  MPI_Helpers::init(argc, argv);

  // oomph_mpi_output.restrict_output_to_single_processor();
#endif

  //--------Store command line arguments-----------------------//
  CommandLineArgs::setup(argc, argv);

  CommandLineArgs::specify_command_line_flag(
    "--restart", &GlobalSimSettings::restart, "Restart");

  CommandLineArgs::specify_command_line_flag(
    "--folder", &GlobalSimSettings::result_folder, "Result Folder");

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
    "--major_axis", &GlobalPhysicalVariables::a, "Ellipse major (x) axis");

  CommandLineArgs::specify_command_line_flag(
    "--minor_axis", &GlobalPhysicalVariables::b, "Ellipse minor (y) axis");

  CommandLineArgs::specify_command_line_flag(
    "--t_max", &GlobalSimSettings::t_max, "Max time");

  CommandLineArgs::specify_command_line_flag(
    "--t_f", &GlobalPhysicalVariables::t_f, "Dryout time");

  CommandLineArgs::specify_command_line_flag(
    "--evaporation_flag",
    &GlobalPhysicalVariables::is_diffusive,
    "Evaporation flag");

  CommandLineArgs::specify_command_line_flag(
    "--initial_height", &GlobalPhysicalVariables::h0, "Initial height");

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
  parameter_file << "a: " << GlobalPhysicalVariables::a << std::endl;
  parameter_file << "b: " << GlobalPhysicalVariables::b << std::endl;
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

  // Create the main problem
  StructuredThinFilmProblem<RefineableQThinFilmBrinkmanElement<2, 3>> problem;

  // Set the initial condition (including initialise dt and impulsive start if
  // we're not restarting a simulation)
  problem.set_initial_condition();

  // Doc initial solution
  problem.doc_solution(true);

  // Want a small time step if we're not restarting, otherwise just take the
  // last dt
  double dt;
  if (GlobalSimSettings::restart)
  {
    dt = problem.time_pt()->time(0) - problem.time_pt()->time(1);
  }
  else
  {
    dt = GlobalSimSettings::dt_comp;
  }

  double time_increment = 0.025;

  // Target time is current time (may not be 0 if we're restarting a
  // simulation)
  // + 0.05
  double t_target = problem.time_pt()->time() + time_increment;

  bool increment_boolean = false;
  bool full_output = false;
  bool success = true;

  double dt_new;

  // If the evaporation is diffusive, will neeed a smaller Pe to begin with
  double Pe_temp = GlobalPhysicalVariables::Pe;
  if (GlobalPhysicalVariables::is_diffusive)
  {
    GlobalPhysicalVariables::Pe = 40.0;
    GlobalSimSettings::adapt_number = 1;
  }
  bool has_been_reset = false;

  //  Timestepping loop:
  unsigned k = 0;
  while (problem.time_pt()->time() < GlobalSimSettings::t_max)
  {
    // If we're nearing the end of the simulation, change the refinement
    // tolerances
    if (problem.time_pt()->time() >= 0.9)
    {
      problem.my_mesh_pt()->max_permitted_error() = 5.0e-5;
      problem.my_mesh_pt()->min_permitted_error() = 1.0e-5;
      GlobalSimSettings::adapt_number = 1;
    }

    // Restore Peclet number to original value
    if ((problem.time_pt()->time() > 0.05) && (!has_been_reset))
    {
      GlobalPhysicalVariables::Pe = Pe_temp;
      has_been_reset = true;
    }

    // If we're running over time (47:00:00), do one final full doc and break
    // out of the loop
    time_t current_time = time(NULL);
    if (current_time - start_time > 169200.0)
    {
      oomph_info
        << "TIME LIMIT HAS BEEN REACHED. BREAKING OUT AND DOING A FINAL "
           "FULL DOC."
        << std::endl;
      problem.doc_solution(true);
      break;
    }

    // Perform solve
    problem.custom_adaptive_unsteady_newton_solve(
      dt, GlobalSimSettings::epsilon_t, dt_new, success);
    if (!success)
    {
      break;
    }

    // The actual timestep taken:
    double dt_taken = problem.time_pt()->time() - problem.time_pt()->time(1);

    // If we are due to increment and the actual timestep was dt, increment
    // t_target by 0.1. Also set the full output boolean to true.
    if (increment_boolean && fabs(dt_taken - dt) < 1.0e-7)
    {
      t_target += time_increment;
      increment_boolean = false;
      full_output = true;
    }

    // If the next trial timestep exceeds the next desired output time,
    // shorten accordingly
    if (problem.time_pt()->time() + dt_new > t_target)
    {
      dt = t_target - problem.time_pt()->time();
      increment_boolean = true;
    }
    else
    {
      dt = dt_new;
    }

    // Doc the solution
    problem.doc_solution(full_output);
    full_output = false;

    // Re-mesh if we're due
    if ((k % GlobalSimSettings::adapt_number == 0) && (k != 0))
    {
      problem.adapt();
    }
    ++k;

    problem.smooth_spurious_nodes();
  }

  // Finalise MPI
#ifdef OOMPH_HAS_MPI

  MPI_Helpers::finalize();

#endif
} // End of main
