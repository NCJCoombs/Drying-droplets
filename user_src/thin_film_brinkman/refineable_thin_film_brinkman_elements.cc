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
#include "refineable_thin_film_brinkman_elements.h"


namespace oomph
{
  //========================================================================
  /// Add element's contribution to the elemental
  /// residual vector and/or Jacobian matrix.
  /// flag=1: compute both
  /// flag=0: compute only residual vector
  //========================================================================
  template<unsigned DIM>
  void RefineableThinFilmBrinkmanEquations<DIM>::
    fill_in_generic_residual_contribution_thin_film_brinkman(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find the index at which the variable is stored
    unsigned u_nodal_index[3];
    for(unsigned a=0; a<3; a++)
    {
     u_nodal_index[a] = this->u_index_thin_film_brinkman(a);
    }

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get the peclet number
    double pe_local = this->peclet();
    
    // Get the dry-out time
    double t_f_local = this->t_f();
    
    // Local variable to determine the ALE stuff
    bool ALE_is_disabled = this->ALE_is_disabled;
    
    // Integers to hold the local equation and unknowns
    int local_eqn = 0, local_unknown = 0;
    
    // Local storage for pointers to hang info objects
    HangInfo *hang_info_pt = 0, *hang_info2_pt = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_thin_film_brinkman(
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
            interpolated_x[j] += nodal_position(l, j) * psi(l);
        }   
          
        // loop over system variables
        for (unsigned a = 0; a < 3; a++)
        {
            // Calculate the value at the nodes
            double u_value = nodal_value(l, u_nodal_index[a]);
            interpolated_u[a] += u_value * psi(l);
            dudt[a] += this->du_dt_thin_film_brinkman(l, a) * psi(l);
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

      // Get source function
      //-------------------
      double evap = 0.0;
      this->get_evap_flux_thin_film_brinkman(time, ipt, interpolated_x, evap);
      
      // Calculate the transition functions and their derivatives
      double transition_term;
      this->get_transition_fct_thin_film_brinkman(ipt, interpolated_u[2], transition_term);
      
      double transition_term_deriv;
      this->get_transition_fct_deriv_thin_film_brinkman(ipt, interpolated_u[2], transition_term_deriv);
      
      // Get mobility
      double mobility;
      this->get_mobility_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], mobility);
      
      // Get the mobility derivatives
      double mobility_deriv_h;
      double mobility_deriv_c;
      this->get_mobility_deriv_h_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], mobility_deriv_h);
      this->get_mobility_deriv_c_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], mobility_deriv_c);
      
      // Get solute mobility
      double solute_mobility;
      this->get_solute_mobility_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], solute_mobility);
      
      // Get the solute mobility derivatives
      double solute_mobility_deriv_h;
      double solute_mobility_deriv_c;
      this->get_solute_mobility_deriv_h_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], solute_mobility_deriv_h);
      this->get_solute_mobility_deriv_c_thin_film_brinkman(ipt, interpolated_u[0], interpolated_u[2], solute_mobility_deriv_c);
      
      // Assemble residuals and Jacobian
      // JACOBIAN IS WRONG, BUT FLAG IS ALWAYS 0 SO DON'T WORRY FOR NOW
      //--------------------------------
      
      bool axisym_flag = this->axisymmetry_flag();
      double r;
      if (axisym_flag)
      {
        r = interpolated_x[0];
      }
      else
      {
        r = 1.0;   
      }

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Local variables to store the number of master nodes and
        // the weight associated with the shape function if the node is hanging
        unsigned n_master = 1;
        double hang_weight = 1.0;
        // Local bool (is the node hanging)
        bool is_node_hanging = this->node_pt(l)->is_hanging();


        // If the node is hanging, get the number of master nodes
        if (is_node_hanging)
        {
          hang_info_pt = this->node_pt(l)->hanging_pt();
          n_master = hang_info_pt->nmaster();
        }
        // Otherwise there is just one master node, the node itself
        else
        {
          n_master = 1;
        }  
        
        // Loop over the number of master nodes
        for (unsigned m = 0; m < n_master; m++)
        {
            ////-------------------CONTRIBUTION FOR h-------------------------------------------------------////
            
            // Get the local equation number and hang_weight
            // If the node is hanging
            if (is_node_hanging)
            {
                // Read out the local equation from the master node
                local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                                u_nodal_index[0]);
                // Read out the weight from the master node
                hang_weight = hang_info_pt->master_weight(m);
            }
            // If the node is not hanging
            else
            {
                // The local equation number comes from the node itself
                local_eqn = this->nodal_local_eqn(l, u_nodal_index[0]);
                // The hang weight is one
                hang_weight = 1.0;
            }

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
                residuals[local_eqn] += test(l) * (pow(transition_term, -1.0) - 1.0) * dudt[0] * r * W * hang_weight;
                residuals[local_eqn] -= test(l) * interpolated_u[1] * r * W * hang_weight;
                
                for (unsigned k = 0; k < DIM; k++)
                {
                  residuals[local_eqn] += dtestdx(l,k) * interpolated_dudx(0,k) * r * W * hang_weight;
                }
                 
              // Calculate the jacobian
              // -----------------------
              if (flag)
              {
                // Local variables to store the number of master nodes
                // and the weights associated with each hanging node
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the nodes for the variables
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local bool (is the node hanging)
                  bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                  // If the node is hanging, get the number of master nodes
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = this->node_pt(l2)->hanging_pt();
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise there is one master node, the node itself
                  else
                  {
                    n_master2 = 1;
                  }
                  
                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    //-----------------------CONTRIBUTION hh---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_index[0]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[0]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                      
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                        jacobian(local_eqn,local_unknown) += test(l) * psi(l2) *
                          node_pt(l2)->time_stepper_pt()->weight(1, 0) * (pow(transition_term, -1.0) - 1.0) * r * W * hang_weight * hang_weight2;
                
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn,local_unknown) += dtestdx(l,k) * dpsidx(l2,k) * r * W * hang_weight * hang_weight2;
                      }    
                    }
                    
                    //-----------------------CONTRIBUTION hw---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_index[1]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[1]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                    
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) -= test(l) * psi(l2) * r * W * hang_weight * hang_weight2;
                    }
                    
                    //-----------------------CONTRIBUTION hc---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                        hang_info2_pt->master_node_pt(m2), u_nodal_index[2]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[2]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }

                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                       jacobian(local_eqn,local_unknown) -= test(l) * pow(transition_term, -2.0) 
                         * transition_term_deriv * psi(l2) * dudt[0] * r * W * hang_weight * hang_weight2;
                    }
                  }
                } // End loop over nodes l2
              } // End of jacobian
            }
                
            ////--------------------CONTRIBUTION FOR w------------------------------------------------------////
            
            // Get the local equation number and hang_weight
            // If the node is hanging
            if (is_node_hanging)
            {
                // Read out the local equation from the master node
                local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                                u_nodal_index[1]);
                // Read out the weight from the master node
                hang_weight = hang_info_pt->master_weight(m);
            }
            // If the node is not hanging
            else
            {
                // The local equation number comes from the node itself
                local_eqn = this->nodal_local_eqn(l, u_nodal_index[1]);
                // The hang weight is one
                hang_weight = 1.0;
            }
            
            /*IF it's not a boundary condition*/
            if (local_eqn >=0)
            {
              // Add evaporation source term and time derivative
              residuals[local_eqn] += (pow(t_f_local, -1.0) * dudt[0] * (1.0 - interpolated_u[2]) 
                                 - pow(t_f_local, -1.0) * interpolated_u[0] * dudt[2] + evap)
                                 * test(l) * r * W * hang_weight;

              // The mesh velocity bit
              if (!ALE_is_disabled)
              {
                for (unsigned k = 0; k < DIM; k++)
                {
                  residuals[local_eqn] -= mesh_velocity[k] *
                                      interpolated_dudx(0,k) * test(l) * r * W * hang_weight;
                }
              }

              // Mobility bit
              for (unsigned k = 0; k < DIM; k++)
              {
                residuals[local_eqn] +=
                 -mobility * interpolated_dudx(1, k) * (1.0 - interpolated_u[2]) * dtestdx(l, k) * r * W * hang_weight;
              }
              
              // Calculate the jacobian
              // -----------------------
              if (flag)
              {
                // Local variables to store the number of master nodes
                // and the weights associated with each hanging node
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the nodes for the variables
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local bool (is the node hanging)
                  bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                  // If the node is hanging, get the number of master nodes
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = this->node_pt(l2)->hanging_pt();
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise there is one master node, the node itself
                  else
                  {
                    n_master2 = 1;
                  }
                  
                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    //-----------------------CONTRIBUTION wh---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[0]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[0]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                
                    //If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        test(l) * psi(l2) * pow(t_f_local, -1.0) *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                        (1.0 - interpolated_u[2]) * r * W * hang_weight * hang_weight2;
                      jacobian(local_eqn, local_unknown) -=
                        test(l) * psi(l2) * pow(t_f_local, -1.0) * dudt[2] * r * W * hang_weight * hang_weight2;
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          dtestdx(l, k) * mobility_deriv_h * psi(l2) *
                        interpolated_dudx(1, k) * (1.0 - interpolated_u[2]) * r * W * hang_weight * hang_weight2;
                      }
                    }
               
                    //-----------------------CONTRIBUTION ww---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[1]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[1]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                  
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn,local_unknown) -= dtestdx(l,k) * mobility * dpsidx(l2,k) * (1.0 - interpolated_u[2]) * r * W * hang_weight * hang_weight2;
                      }
                    }
                  
                    //-----------------------CONTRIBUTION wc---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[2]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[2]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }

                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        test(l) * pow(t_f_local, -1.0) * dudt[0] * psi(l2) * r * W * hang_weight * hang_weight2;
                      jacobian(local_eqn, local_unknown) -=
                        test(l) * pow(t_f_local, -1.0) * interpolated_u[0] * psi(l2) *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0) * r * W * hang_weight * hang_weight2;

                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                            dtestdx(l, k) * mobility_deriv_c * psi(l2) *
                            interpolated_dudx(1, k) * (1.0 - interpolated_u[2]) * r * W * hang_weight * hang_weight2;
                        jacobian(local_eqn, local_unknown) +=
                            dtestdx(l, k) * mobility * interpolated_dudx(1, k) *
                            psi(l2) * r * W * hang_weight * hang_weight2;
                      }
                    }
                  } // End of master node loop m2
                } // End of node loop l2
              } // End of Jacobian
            }

            ////-------------------CONTRIBUTION FOR c-------------------------------------------------------////
            
            // Get the local equation number and hang_weight
            // If the node is hanging
            if (is_node_hanging)
            {
                // Read out the local equation from the master node
                local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                                u_nodal_index[2]);
                // Read out the weight from the master node
                hang_weight = hang_info_pt->master_weight(m);
            }
            // If the node is not hanging
            else
            {
                // The local equation number comes from the node itself
                local_eqn = this->nodal_local_eqn(l, u_nodal_index[2]);
                // The hang weight is one
                hang_weight = 1.0;
            }
            
            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Add time derivatives
              residuals[local_eqn] += pow(t_f_local, -1.0) * (dudt[0]*interpolated_u[2] + dudt[2]*interpolated_u[0])
                                      * test(l) * r * W * hang_weight;

              // The mesh velocity bit
              if (!ALE_is_disabled)
              {
                for (unsigned k = 0; k < DIM; k++)
                {
                  residuals[local_eqn] -= mesh_velocity[k] *
                                      (interpolated_dudx(0,k)*interpolated_u[2] + interpolated_dudx(2,k)*interpolated_u[0])
                                      * test(l) * r * W * hang_weight;
                }
              }

              // Mobility bit
              for (unsigned k = 0; k < DIM; k++)
              {
                residuals[local_eqn] +=
                -solute_mobility * interpolated_dudx(1, k) * interpolated_u[2] * dtestdx(l, k) * r * W * hang_weight;
              }
          
              // The diffusion bit
              for (unsigned k = 0; k < DIM; k++)
              {
                residuals[local_eqn] +=
                  interpolated_u[0] * dtestdx(l,k) * interpolated_dudx(2,k) * pow(pe_local,-1.0) * r * W * hang_weight;
              }
              
              // Calculate the jacobian
              // -----------------------
              if (flag)
              {
                // Local variables to store the number of master nodes
                // and the weights associated with each hanging node
                unsigned n_master2 = 1;
                double hang_weight2 = 1.0;
                // Loop over the nodes for the variables
                for (unsigned l2 = 0; l2 < n_node; l2++)
                {
                  // Local bool (is the node hanging)
                  bool is_node2_hanging = this->node_pt(l2)->is_hanging();
                  // If the node is hanging, get the number of master nodes
                  if (is_node2_hanging)
                  {
                    hang_info2_pt = this->node_pt(l2)->hanging_pt();
                    n_master2 = hang_info2_pt->nmaster();
                  }
                  // Otherwise there is one master node, the node itself
                  else
                  {
                    n_master2 = 1;
                  }
                  
                  // Loop over the master nodes
                  for (unsigned m2 = 0; m2 < n_master2; m2++)
                  {
                    //-----------------------CONTRIBUTION ch---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[0]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[0]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                
                    //If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        test(l) * psi(l2) * pow(t_f_local, -1.0) *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0) *
                        interpolated_u[2] * r * W * hang_weight * hang_weight2;

                      jacobian(local_eqn, local_unknown) +=
                        test(l) * psi(l2) * pow(t_f_local, -1.0) * dudt[2] * r * W * hang_weight * hang_weight2;

                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                        dtestdx(l, k) *
                        (solute_mobility_deriv_h * psi(l2) *
                        interpolated_dudx(1, k) * interpolated_u[2] -
                        pow(pe_local, -1.0) * psi(l2) * interpolated_dudx(2, k)) *
                        r * W * hang_weight * hang_weight2;
                      }
                    }
               
                    //-----------------------CONTRIBUTION cw---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[1]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[1]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }
                  
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      for (unsigned k = 0; k < DIM; k++)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          dtestdx(l, k) * solute_mobility * dpsidx(l2, k) *
                          interpolated_u[2] * r * W * hang_weight * hang_weight2;
                      }
                    }
                  
                    //-----------------------CONTRIBUTION cc---------------------//
                    // Get the local unknown and weight
                    // If the node is hanging
                    if (is_node2_hanging)
                    {
                      // Read out the local unknown from the master node
                      local_unknown = this->local_hang_eqn(
                      hang_info2_pt->master_node_pt(m2), u_nodal_index[2]);
                      // Read out the hanging weight from the master node
                      hang_weight2 = hang_info2_pt->master_weight(m2);
                    }
                    // If the node is not hanging
                    else
                    {
                      // The local unknown number comes from the node
                      local_unknown =
                        this->nodal_local_eqn(l2, u_nodal_index[2]);
                      // The hang weight is one
                      hang_weight2 = 1.0;
                    }

                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        test(l) * pow(t_f_local, -1.0) *
                        (dudt[0] * psi(l2) +
                        interpolated_u[0] * psi(l2) *
                        node_pt(l2)->time_stepper_pt()->weight(1, 0)) *
                        r * W * hang_weight * hang_weight2;

                      for (unsigned k = 0; k < DIM; k++)
                      {
                       jacobian(local_eqn, local_unknown) -=
                        dtestdx(l, k) *
                        (solute_mobility_deriv_c * interpolated_u[2] +
                        solute_mobility) *
                        psi(l2) * interpolated_dudx(1, k) * r * W * hang_weight * hang_weight2;
                        
                       jacobian(local_eqn, local_unknown) +=
                        dtestdx(l, k) * pow(pe_local, -1.0) * interpolated_u[0] *
                        dpsidx(l2, k) * r * W * hang_weight * hang_weight2;
                      }
                    }
                  } // End of master node loop m2
                } // End of node loop l2
              } // End of Jacobian
            }
        } // End of loop over master nodes
      } // End of loop over test functions
    } // End of loop over integration points
  }

  //====================================================================
  // Force build of templates
  //====================================================================
  template class RefineableQThinFilmBrinkmanElement<1, 2>;
  template class RefineableQThinFilmBrinkmanElement<1, 3>;
  template class RefineableQThinFilmBrinkmanElement<1, 4>;
  
  template class RefineableQThinFilmBrinkmanElement<2, 2>;
  template class RefineableQThinFilmBrinkmanElement<2, 3>;
  template class RefineableQThinFilmBrinkmanElement<2, 4>;

} // namespace oomph
