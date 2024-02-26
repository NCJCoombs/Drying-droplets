// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Non-inline functions for UnsteadyHeat elements
#include "thin_film_brinkman_elements.h"


namespace oomph
{
  /// 1D/2D Thin film Brinkman elements

  /// Default value for the Peclet number
  template<unsigned DIM>
  double ThinFilmBrinkmanEquations<DIM>::Default_peclet_number = 1.0;

  /// Default value for the dry-out time
  template<unsigned DIM>
  double ThinFilmBrinkmanEquations<DIM>::Default_dryout_time = 1.0;

  /// Default axisymmetry flag
  template<unsigned DIM>
  bool ThinFilmBrinkmanEquations<DIM>::Default_axisymmetry_flag = false;


  //======================================================================
  // Set the data for the number of Variables at each node
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  const unsigned QThinFilmBrinkmanElement<DIM, NNODE_1D>::Initial_Nvalue = 3;

  //======================================================================
  /// Compute element residual Vector and/or element Jacobian matrix
  ///
  /// flag=1: compute both
  /// flag=0: compute only residual Vector
  ///
  /// Pure version without hanging nodes
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::
    fill_in_generic_residual_contribution_thin_film_brinkman(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find the index at which the variable is stored
    unsigned u_nodal_index[3];
    for (unsigned a = 0; a < 3; a++)
    {
      u_nodal_index[a] = u_index_thin_film_brinkman(a);
    }

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get the peclet number
    double pe_local = peclet();

    // Get the dry-out time
    double t_f_local = t_f();

    // Integers to hold the local equation and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = dshape_and_dtest_eulerian_at_knot_thin_film_brinkman(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate memory for local quantities and initialise to zero
      Vector<double> interpolated_u(3, 0.0);
      Vector<double> dudt(3, 0.0);
      Vector<double> interpolated_x(DIM, 0.0);
      DenseMatrix<double> interpolated_dudx(3, DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);

      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
        }

        // loop over system variables
        for (unsigned a = 0; a < 3; a++)
        {
          // Calculate the value at the nodes
          double u_value = raw_nodal_value(l, u_nodal_index[a]);
          interpolated_u[a] += u_value * psi(l);
          dudt[a] += du_dt_thin_film_brinkman(l, a) * psi(l);
          // Loop over directions
          for (unsigned j = 0; j < DIM; j++)
          {
            interpolated_dudx(a, j) += u_value * dpsidx(l, j);
          }
        }
      }

      // Mesh velocity?
      if (!ALE_is_disabled)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            mesh_velocity[j] += raw_dnodal_position_dt(l, j) * psi(l);
          }
        }
      }

      // Get evaporation flux. Either use the vector pointer if it has
      // been assigned, or use the access function
      //---------------------------------------------------------------
      double evap = 0.0;
      if (Evap_at_int_points_pt)
      {
        Vector<double> evap_temp(n_intpt, 0.0);
        evap_temp = evap_at_int_points();
        evap = evap_temp[ipt];
      }
      else
      {
        get_evap_flux_thin_film_brinkman(time, ipt, interpolated_x, evap);
      }

      // Calculate the transition functions and their derivatives
      double transition_term;
      get_transition_fct_thin_film_brinkman(
        ipt, interpolated_u[2], transition_term);

      double transition_term_deriv;
      get_transition_fct_deriv_thin_film_brinkman(
        ipt, interpolated_u[2], transition_term_deriv);

      // Get mobility
      double mobility;
      get_mobility_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], mobility);

      // Get the mobility derivatives
      double mobility_deriv_h;
      double mobility_deriv_c;
      get_mobility_deriv_h_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], mobility_deriv_h);
      get_mobility_deriv_c_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], mobility_deriv_c);

      // Get solute mobility
      double solute_mobility;
      get_solute_mobility_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], solute_mobility);

      // Get the solute mobility derivatives
      double solute_mobility_deriv_h;
      double solute_mobility_deriv_c;
      get_solute_mobility_deriv_h_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], solute_mobility_deriv_h);
      get_solute_mobility_deriv_c_thin_film_brinkman(
        ipt, interpolated_u[0], interpolated_u[2], solute_mobility_deriv_c);

      // Get the axisymmetry flag
      bool axisym_flag = axisymmetry_flag();
      double r;
      if (axisym_flag)
      {
        r = interpolated_x[0];
      }
      else
      {
        r = 1.0;
      }

      // Assemble residuals and Jacobian

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        ////--------------------CONTRIBUTION FOR
        /// h------------------------------------------------------////

        local_eqn = nodal_local_eqn(l, u_nodal_index[0]);

        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          residuals[local_eqn] +=
            test(l) * (pow(transition_term, -1.0) - 1.0) * dudt[0] * r * W;
          residuals[local_eqn] -= test(l) * interpolated_u[1] * r * W;

          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] +=
              dtestdx(l, k) * interpolated_dudx(0, k) * r * W;
          }

          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              //-----------------------CONTRIBUTION hh---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  test(l) * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                  (pow(transition_term, -1.0) - 1.0) * r * W;

                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) +=
                    dtestdx(l, k) * dpsidx(l2, k) * r * W;
                }
              }

              //-----------------------CONTRIBUTION hw---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -= test(l) * psi(l2) * r * W;
              }

              //-----------------------CONTRIBUTION hc---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -=
                  test(l) * pow(transition_term, -2.0) * transition_term_deriv *
                  psi(l2) * dudt[0] * r * W;
              }
            }
          } // End of jacobian
        }

        ////-------------------CONTRIBUTION FOR
        /// w-------------------------------------------------------////

        local_eqn = nodal_local_eqn(l, u_nodal_index[1]);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add evaporation source term and time derivative
          residuals[local_eqn] +=
            (pow(t_f_local, -1.0) * dudt[0] * (1.0 - interpolated_u[2]) -
             pow(t_f_local, -1.0) * interpolated_u[0] * dudt[2] + evap) *
            test(l) * r * W;

          // The mesh velocity bit
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -=
                mesh_velocity[k] * interpolated_dudx(0, k) * test(l) * r * W;
            }
          }

          // Mobility bit
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] += -mobility * interpolated_dudx(1, k) *
                                    (1.0 - interpolated_u[2]) * dtestdx(l, k) *
                                    r * W;
          }

          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              //-----------------------CONTRIBUTION wh---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  test(l) * psi(l2) * pow(t_f_local, -1.0) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                  (1.0 - interpolated_u[2]) * r * W;
                jacobian(local_eqn, local_unknown) -=
                  test(l) * psi(l2) * pow(t_f_local, -1.0) * dudt[2] * r * W;

                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) * mobility_deriv_h * psi(l2) *
                    interpolated_dudx(1, k) * (1.0 - interpolated_u[2]) * r * W;
                }
              }

              //-----------------------CONTRIBUTION ww---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) * mobility * dpsidx(l2, k) *
                    (1.0 - interpolated_u[2]) * r * W;
                }
              }

              //-----------------------CONTRIBUTION wc---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) -=
                  test(l) * pow(t_f_local, -1.0) * dudt[0] * psi(l2) * r * W;
                jacobian(local_eqn, local_unknown) -=
                  test(l) * pow(t_f_local, -1.0) * interpolated_u[0] * psi(l2) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) * r * W;

                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) * mobility_deriv_c * psi(l2) *
                    interpolated_dudx(1, k) * (1.0 - interpolated_u[2]) * r * W;
                  jacobian(local_eqn, local_unknown) +=
                    dtestdx(l, k) * mobility * interpolated_dudx(1, k) *
                    psi(l2) * r * W;
                }
              }
            }
          }
        }

        ////-------------------CONTRIBUTION FOR
        /// c-------------------------------------------------------////

        local_eqn = nodal_local_eqn(l, u_nodal_index[2]);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add time derivatives
          residuals[local_eqn] +=
            pow(t_f_local, -1.0) *
            (dudt[0] * interpolated_u[2] + dudt[2] * interpolated_u[0]) *
            test(l) * r * W;

          // The mesh velocity bit
          if (!ALE_is_disabled)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -=
                mesh_velocity[k] *
                (interpolated_dudx(0, k) * interpolated_u[2] +
                 interpolated_dudx(2, k) * interpolated_u[0]) *
                test(l) * r * W;
            }
          }

          // Mobility bit
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] += -solute_mobility * interpolated_dudx(1, k) *
                                    interpolated_u[2] * dtestdx(l, k) * r * W;
          }

          // The diffusion bit
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] += interpolated_u[0] * dtestdx(l, k) *
                                    interpolated_dudx(2, k) *
                                    pow(pe_local, -1.0) * r * W;
          }

          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              //-----------------------CONTRIBUTION ch---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[0]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  test(l) * psi(l2) * pow(t_f_local, -1.0) *
                  node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                  interpolated_u[2] * r * W;

                jacobian(local_eqn, local_unknown) +=
                  test(l) * psi(l2) * pow(t_f_local, -1.0) * dudt[2] * r * W;

                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) *
                    (solute_mobility_deriv_h * psi(l2) *
                       interpolated_dudx(1, k) * interpolated_u[2] -
                     pow(pe_local, -1.0) * psi(l2) * interpolated_dudx(2, k)) *
                    r * W;
                }
              }

              //-----------------------CONTRIBUTION cw---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[1]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) * solute_mobility * dpsidx(l2, k) *
                    interpolated_u[2] * r * W;
                }
              }

              //-----------------------CONTRIBUTION cc---------------------//

              // Set the number of the unknown
              local_unknown = nodal_local_eqn(l2, u_nodal_index[2]);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                jacobian(local_eqn, local_unknown) +=
                  test(l) * pow(t_f_local, -1.0) *
                  (dudt[0] * psi(l2) +
                   interpolated_u[0] * psi(l2) *
                     node_pt(l2)->time_stepper_pt()->weight(1, 0)) *
                  r * W;

                for (unsigned k = 0; k < DIM; k++)
                {
                  jacobian(local_eqn, local_unknown) -=
                    dtestdx(l, k) *
                    (solute_mobility_deriv_c * interpolated_u[2] +
                     solute_mobility) *
                    psi(l2) * interpolated_dudx(1, k) * r * W;
                  jacobian(local_eqn, local_unknown) +=
                    dtestdx(l, k) * pow(pe_local, -1.0) * interpolated_u[0] *
                    dpsidx(l2, k) * r * W;
                }
              }
            }
          }
        }
      } // End of loop over test functions
    } // End of loop over integration points
  }


  //======================================================================
  /// Compute norm of fe solution
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::compute_norm(double& norm)
  {
    // Initialise
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get FE function value: just use the height
      double u = interpolated_u_thin_film_brinkman(s, 0);

      // Add to  norm
      norm += u * u * W;
    }
  }

  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  template<unsigned DIM>
  unsigned ThinFilmBrinkmanEquations<DIM>::self_test()
  {
    bool passed = true;

    // Check lower-level stuff
    if (FiniteElement::self_test() != 0)
    {
      passed = false;
    }

    // Return verdict
    if (passed)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }


  //======================================================================
  /// Output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::output(std::ostream& outfile,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    //     Tecplot header info
    //     outfile << tecplot_zone_string(nplot);
    
    // Find out how many nodes there are
    unsigned n_node = nnode();
    
    // Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node, DIM);        

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);
      
      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psif, dpsifdx);
      
      // Allocate for pressure derivates
      Vector<double> interpolated_dpdx(DIM, 0.0);

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions for velocity gradients
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_dpdx[j] +=
              nodal_value(l, 1) * dpsifdx(l, j);
        }
      }
      
      // Calculate the (solvent) velocity field:
      unsigned dummy_ipt = 0;
      double h_local = interpolated_u_thin_film_brinkman(s, 0);
      double phi_local = interpolated_u_thin_film_brinkman(s, 2);
      Vector<double> vel_field(DIM, 0.0);
      double mobility;
      get_mobility_thin_film_brinkman(
        dummy_ipt, h_local, phi_local, mobility);
      
      for (unsigned j = 0; j < DIM; j++)
      {
        vel_field[j] = mobility * interpolated_dpdx[j] / h_local;
      }

      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << interpolated_x(s, i) << " ";
      }
      for (unsigned i = 0; i < 3; i++)
      {
        outfile << interpolated_u_thin_film_brinkman(s, i) << " ";
       // outfile << interpolated_du_dt_thin_film_brinkman(s, i) << " ";
      }
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << vel_field[i] << " ";
      }
      outfile << std::endl;
    }

    //     Write tecplot footer (e.g. FE connectivity lists)
    //     write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::output(FILE* file_pt,
                                              const unsigned& nplot)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Tecplot header info
    fprintf(file_pt, "%s", tecplot_zone_string(nplot).c_str());

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      for (unsigned i = 0; i < DIM; i++)
      {
        fprintf(file_pt, "%g ", interpolated_x(s, i));
      }
      for (unsigned a = 0; a < 3; a++)
      {
        fprintf(file_pt, "%g \n", interpolated_u_thin_film_brinkman(s, a));
      }
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(file_pt, nplot);
  }


  //======================================================================
  /// Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Output exact solution at time t
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    const double& time,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)

  {
    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


    // Tecplot header info
    outfile << tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over plot points
    unsigned num_plot_points = nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);

    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot header info
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      double u_fe = interpolated_u_thin_film_brinkman(s, 0);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //======================================================================
  /// Validate against exact solution at time t.
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  template<unsigned DIM>
  void ThinFilmBrinkmanEquations<DIM>::compute_error(
    std::ostream& outfile,
    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
    const double& time,
    double& error,
    double& norm)

  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Vector of local coordinates
    Vector<double> s(DIM);

    // Vector for coordintes
    Vector<double> x(DIM);


    // Find out how many nodes there are in the element
    unsigned n_node = nnode();

    Shape psi(n_node);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Tecplot
    outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(1);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++)
      {
        s[i] = integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J = J_eulerian(s);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      double u_fe = interpolated_u_thin_film_brinkman(s, 0);

      // Get exact solution at this point
      (*exact_soln_pt)(time, x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < DIM; i++)
      {
        outfile << x[i] << " ";
      }
      outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

      // Add to error and norm
      norm += exact_soln[0] * exact_soln[0] * W;
      error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
    }
  }


  //====================================================================
  // Force build of templates (1 and 2 dimensional problems)
  //====================================================================
  template class QThinFilmBrinkmanElement<1, 2>;
  template class QThinFilmBrinkmanElement<1, 3>;
  template class QThinFilmBrinkmanElement<1, 4>;

  template class QThinFilmBrinkmanElement<2, 2>;
  template class QThinFilmBrinkmanElement<2, 3>;
  template class QThinFilmBrinkmanElement<2, 4>;

} // namespace oomph
