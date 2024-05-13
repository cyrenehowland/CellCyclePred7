/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
double max_life_cycle_p = 1.0;
double min_life_cycle_p = 0.0;
const double EPSILON = 1e-6;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~ Setup Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




// =============================================================================
// Environment Setup Functions
// =============================================================================

void setup_microenvironment( void )
{
    // set domain parameters
    
    // put any custom code to set non-homogeneous initial conditions or
    // extra Dirichlet nodes here.
    
    // initialize BioFVM
    
    initialize_microenvironment();
    
//    int idx_food = microenvironment.find_density_index( "food" );
//    
//    // Gaussian distribution parameters
//    double mean_x_A = microenvironment.mesh.bounding_box[3] / 2;  // center x
//    double mean_y_A = microenvironment.mesh.bounding_box[4] / 2;  // center y
//    double sigma_A = 200.0; // standard deviation in microns
//    double amplitude_A = 100.0; // Amplitude of the Gaussian distribution
// 
//    for( int n=0; n < microenvironment.number_of_voxels(); n++ )
//    {
//        auto& voxel = microenvironment.voxels(n);
//                double x_A = voxel.center[0];
//                double y_A = voxel.center[1];
//                
//                // Calculate the Gaussian function
//                double value_A = amplitude_A * exp( - (pow(x_A - mean_x_A, 2) + pow(y_A - mean_y_A, 2)) / (2 * pow(sigma_A, 2)) );
//                microenvironment(n)[idx_food] = value_A;  // assuming the substrate index is 0
//
//    }
//   
//    store_initial_concentrations();
//    return;
}





// =============================================================================
// Cell Type Definitions, Initializations, Setup and Tissue Configuration
// =============================================================================

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
    
    // Make sure we're ready for 2D
    cell_defaults.functions.set_orientation = up_orientation;
    cell_defaults.phenotype.geometry.polarity = 1.0;
    cell_defaults.phenotype.motility.restrict_to_2D = true;
    
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ~~~~~~~~~~~~~ (Prey) Cell Function Assignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cell_Definition* pPreyDef = find_cell_definition( "Cell" );
    
//    pPreyDef -> functions.update_velocity = prey_update_cell_velocity;
//    pPreyDef -> functions.update_migration_bias = prey_migration_bias_function;
    
    pPreyDef->functions.update_phenotype = prey_phenotype_function;
    pPreyDef->phenotype.mechanics.attachment_elastic_constant = parameters.doubles("attachment_elastic_constant");
    pPreyDef->functions.contact_function = standard_elastic_contact_function;
    
    // Cell Cycle:
    
    Cycle_Model& prey_cycle_model = pPreyDef -> phenotype.cycle.model(); // Access cycle model for modification

    // Don't let division occur automatically after any phase
    for (int i = 0; i < prey_cycle_model.phases.size(); ++i)
    {
        prey_cycle_model.phases[i].division_at_phase_exit = false; // reset division flag for all phases
    }
//    prey_cycle_model.phases[1].division_at_phase_exit = true;

    prey_cycle_model.phases[0].entry_function = phase_0_entry_function;
    prey_cycle_model.phases[1].entry_function = phase_1_entry_function;
    prey_cycle_model.phases[2].entry_function = phase_2_entry_function;
    prey_cycle_model.phases[3].entry_function = NULL;
    
    prey_cycle_model.phase_link(0,1).exit_function = phase_0_exit_function;
    prey_cycle_model.phase_link(1,2).exit_function = phase_1_exit_function;
    prey_cycle_model.phase_link(2,3).exit_function = phase_2_exit_function;
    prey_cycle_model.phase_link(3,0).exit_function = phase_3_exit_function;
    
    prey_cycle_model.phase_link(0,1).arrest_function = G0_arrest_function; // Assign an arrest function
//    prey_cycle_model.phase_link(0,1).arrest_function = NULL;
    prey_cycle_model.phase_link(1,2).arrest_function = S_arrest_function;
    prey_cycle_model.phase_link(2,3).arrest_function = NULL;
    prey_cycle_model.phase_link(3,0).arrest_function = NULL;
        
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ~~~~~~~~~~~~~ Predator Cell Function Assignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cell_Definition* pPredDef = find_cell_definition( "Predator" );
    
    Cycle_Model& pred_cycle_model = pPredDef -> phenotype.cycle.model(); // Access cycle model for modification
    
    pred_cycle_model.phase_link(0,1).arrest_function = pred_repo_arrest_function; // Assign an arrest function

    
    pPredDef->functions.update_phenotype = pred_phenotype_function;
    
	display_cell_definitions( std::cout );
	
	return; 
}


//void initialize_random_orientation(Cell* pCell)
//{
//    // Generate random angles in radians
//    double theta = UniformRandom() * 2.0 * M_PI; // Random angle from 0 to 2*pi
//   
//    // Convert polar coordinates to Cartesian coordinates for the unit vector
//    double x = cos(theta);
//    double y = sin(theta);
//
//    // Assign to cell's orientation
//    pCell->state.orientation[0] = x;
//    pCell->state.orientation[1] = y;
//    pCell->state.orientation[2] = 0;
//}


void prey_initialize_properties(Cell* cell){
    cell -> custom_data["initial_volume"] = cell -> phenotype.volume.total;
    cell -> phenotype.motility.is_motile = true;
    cell -> custom_data["energy"] = 25;
    cell -> custom_data["stress"] = 0;

    
    
    // Set cell cycle transition rates:

    // Randomizing initially life_cycle_p for each cell
    cell->custom_data["life_cycle_p"] = UniformRandom() * (max_life_cycle_p - min_life_cycle_p) + min_life_cycle_p;

    // Retrieve custom data from the specific cell
    double life_history_t = cell->custom_data["life_history_t"];
    double life_cycle_p = cell->custom_data["life_cycle_p"];
    
    // Set transition rates for phase 0 to 1 and phase 2 to 3
    double phase_0_duration = (1 - life_cycle_p) * life_history_t;
    double phase_2_duration = life_cycle_p * life_history_t;
    
    // Use EPSILON to ensure durations are not zero
    phase_0_duration = std::max(phase_0_duration, EPSILON);
    phase_2_duration = std::max(phase_2_duration, EPSILON);

    
    // Modify the cell cycle of the specific cell
    cell->phenotype.cycle.data.transition_rate(0, 1) = 1.0 / phase_0_duration;
    cell->phenotype.cycle.data.transition_rate(2, 3) = 1.0 / phase_2_duration;
    
    // Initialize random orientation
//    initialize_random_orientation(cell);

}


void pred_initialize_properties(Cell* cell){
    cell -> custom_data["energy"] = 10;
    cell -> phenotype.motility.is_motile = true;

    
}


void setup_tissue( void )
{
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];
    
    if( default_microenvironment_options.simulate_2D == true )
    {
        Zmin = 0.0;
        Zmax = 0.0;
    }
    
    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;
    
    // create some of each type of cell
    
    Cell* pC;
    
    for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    {
        Cell_Definition* pCD = cell_definitions_by_index[k];
        std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
        for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
        {
            std::vector<double> position = {0,0,0};
            position[0] = Xmin + UniformRandom()*Xrange;
            position[1] = Ymin + UniformRandom()*Yrange;
            position[2] = Zmin + UniformRandom()*Zrange;
            
            pC = create_cell( *pCD );
            pC->assign_position( position );
        }
    }
    std::cout << std::endl;
    
    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();
    
    // Initialize cell Properties:
    for( auto& cell : *all_cells )
    {
        // Check if the current cell is a default cell
        if (cell->type == find_cell_definition("Cell")->type){
            
            // Initilaize cell properties:
            prey_initialize_properties(cell);
        }
        
        if (cell->type == find_cell_definition("Predator")->type){
            
            // Initilaize cell properties:
            pred_initialize_properties(cell);
        }
    }
    
    return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~ Update Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




// =============================================================================
// Environment Update Function
// =============================================================================

//void update_microenvironment(double dt)
//{
//    static int food_idx = microenvironment.find_density_index("food");
//    double regeneration_rate = 0.1; // set this to your desired rate
//
//    #pragma omp parallel for
//    for(int i = 0; i < microenvironment.mesh.voxels.size(); i++)
//    {
//        auto& voxel = microenvironment.mesh.voxels[i];
//        double updated_concentration = microenvironment.density_vector(i)[food_idx] + regeneration_rate * dt;
//        // Ensure the concentration does not exceed the initial value
//        if(updated_concentration > initial_food_concentration[i])
//        {
//            microenvironment.density_vector(i)[food_idx] = initial_food_concentration[i];
//        }
//        else
//        {
//            microenvironment.density_vector(i)[food_idx] = updated_concentration;
//        }
//    }
//    
//    
//}



// =============================================================================
// Core PhysiCell Cell Update Functions
// =============================================================================

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 


// ****************************************************************************
// ****************************************************************************
// PREY FUNCTONS
// ****************************************************************************
// ****************************************************************************

void prey_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    // If you are attached to others, label 'inedible'
    if (pCell->state.number_of_attached_cells() > 0){
        pCell -> custom_data["edible"] = false;
    }
    else {pCell -> custom_data["edible"] = true;}
    
    
//    std::cout << "Cell (ID: " << pCell->ID << " ) has number of attachments =  " << pCell -> state.number_of_attached_cells() << std::endl;
    
    // Make sure growth only happens in the first and third cycle phases
    int current_phase_index = pCell->phenotype.cycle.current_phase_index();
    if ( current_phase_index == 0 || current_phase_index == 2){
        prey_growth_and_metabolism(pCell, phenotype, dt);
    }
    
    // Stress acounting based on enviornment (low food or predator presence
    prey_stress_function(pCell, phenotype, dt);
    
    // Death based on starvation
    if( pCell->custom_data["energy"] < 1.0 )
    {
        pCell->lyse_cell();

    }

    return; }

// =============================================================================
//  Prey Motility and Motion Update Functions
// =============================================================================

//void custom_update_motility_vector(Cell* pCell, Phenotype& phenotype , double dt)
//{
//    double orientation_bias;
//    
//    if( phenotype.motility.is_motile == false )
//    {
//        phenotype.motility.motility_vector.assign( 3, 0.0 );
//        return;
//    }
//    
//    if( UniformRandom() < dt / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt )
//    {
//        std::vector<double> randvec(3,0.0);
//        if( phenotype.motility.restrict_to_2D == true )
//        { randvec = UniformOnUnitCircle(); }
//        else
//        { randvec = UniformOnUnitSphere(); }
//        
//        // if the update_bias_vector function is set, use it
//        if( pCell -> functions.update_migration_bias )
//        {
//            pCell -> functions.update_migration_bias( pCell ,phenotype,dt );
//        }
//        
//        // Calculate new bias direction taking into account orientation and chemical gradient:
//        
//        // Fetch orientation vector
//        std::vector<double> orient_vector = pCell -> state.orientation;
//        
//        // Fetch chemotaxis vector
//        static int n_index = microenvironment.find_density_index( "food" );
//        std::vector<double> chemo_grad = pCell->nearest_gradient(n_index);
//        
//        normalize( &chemo_grad ); // Normalize
//        
//        
//        if ( norm(chemo_grad) > 1e-3 )
//            {
//                phenotype.motility.migration_bias = std::min(norm(chemo_grad) + 0.6, 0.99 );
//                
//                // Take dot product of both normalized vectors
//                double dot_prod = dot_product(orient_vector, chemo_grad);
//                
////                std::cout << "Dot product of vectors for Cell (ID: " << pCell-> ID << " ) : " << dot_prod << std::endl;
//                
//                // Use dot product to determine new direction
//                if (dot_prod > 0.5) // Vectors are aligned
//                {
//                    phenotype.motility.persistence_time = 2.0 + 4*dot_prod;
//                }
//                else // Vectors are not aligned
//                {
//                    // If the cell is facing away from the detected gradient direction, update orientation more frequently
//                    phenotype.motility.persistence_time  = 2.0;
//                }
//            }
//            else {
//                phenotype.motility.persistence_time = 5.0;
//                phenotype.motility.migration_speed = 5.0;
//                phenotype.motility.migration_bias = 0.6;
//            }
//        
//        
//        
//        
//        // Set orientation bias:
//        // (how much the cell corrects orientation towards chemical gradient direction after each run)
//        // How much the cells movement is determined by its orientatioin vs the detected chemical gradient direction
//        orientation_bias = 0.5;
//        
//        // Update migration bias direction:
//        phenotype.motility.migration_bias_direction = orientation_bias*chemo_grad; //bias_vector = orientation_bias*chem_grad_direction
//
//        double one_minus_orientation_bias = 1.0 - orientation_bias;
//        axpy(&(phenotype.motility.migration_bias_direction), one_minus_orientation_bias, orient_vector); //bias_vector = orientation_bias*chem_grad_direction + (1-orientation_bias)*orientation_direction
//        normalize( &(phenotype.motility.migration_bias_direction) );
//
//
//        
//        // Update motility vector:
//        phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motility = bias_vector
//        phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector
//        
//        double one_minus_bias = 1.0 - (phenotype.motility.migration_bias);
//        
//        axpy( &phenotype.motility.motility_vector, one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
//        normalize( &phenotype.motility.motility_vector );
//        
//        // Update cell's new orientation to new motility vector direction (normalized)
////        pCell -> state.orientation = phenotype.motility.motility_vector;
//        // Scale motility vector by migration speed
//        phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
//        
//        
////        phenotype.motility.motility_vector[2] = 0.0; // Enforce no movement in z-direction
////        pCell->state.orientation[2] = 0.0; // Enforce no orientation component in z-direction
//    }
//    return;
//}


//void custom_migration_bias_function( Cell* pCell, Phenotype& phenotype , double dt ){
//    
//    // Fetch chemotaxis vector
//    static int food_idx = microenvironment.find_density_index( "food" );
//    std::vector<double> chemo_grad = pCell->nearest_gradient(food_idx);
//    
//    
//    std::cout << "norm of chemo gradient for Cell (ID: " << pCell-> ID << ") is = " << norm(chemo_grad) << std::endl;
//
//        if (norm(chemo_grad) > 0.1){
//            phenotype.motility.migration_bias = std::min(norm(chemo_grad) + 0.6, 0.99 );
//
//        }
//        else {
//            phenotype.motility.migration_bias = 0.6;
//
//        }


// =============================================================================
// Prey Life Cycle and Environmental Interaction Management
// =============================================================================

void prey_growth_and_metabolism(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Access index for food
    int food_index = microenvironment.find_density_index("food");
    
    // sample the microenvironment at the cell’s locaiton
    double food_nearby = pCell->nearest_density_vector()[food_index];
    
    // Define the threshold and calculate food_eaten
    double threshold = 10.0;  // Example threshold
    double food_eaten;
    if (food_nearby < threshold)
        food_eaten = food_nearby;  // Eat all available food if below threshold
    else
        food_eaten = 10.0;  // Eat a fixed amount otherwise
       
    
    // Update Energy:
    // Cell gains energy equivalent to food consumption
    pCell->custom_data["energy"] += food_eaten; // Cell gains energy equivalent to food consumption
    // Update energy based on metabolism
    static double metabolic_rate = pCell -> custom_data["metabolic_rate"];
    pCell->custom_data["energy"] /= (1.0 + dt*metabolic_rate);
    
    // Update the microenvironment:
    // reduce food at the cell's location (it "ate" it)
    double new_food_level = std::max(0.0, food_nearby - food_eaten); // Ensure not negative
    int voxel_index = pCell->get_current_voxel_index();
    microenvironment(voxel_index)[food_index] = new_food_level;
    

    
    // Update Volume:
    // Constants for the power law model
    const double b = 3/4 ; // Define the exponent b
    const double k_base = 10;  // Base value for k
    double k = k_base * food_eaten;
    
    // Get the current volume
    double current_volume = pCell->phenotype.volume.total;
    // Calculate the rate of change of volume
    double dVdt = k * pow(current_volume, b);
    // Update the volume based on the rate of change and the time step
    double new_volume = current_volume + dVdt * dt;
    // Set the new volume
    pCell->set_total_volume(new_volume);
    
    //    std::cout << "Cell ID ("<< pCell->ID <<" ) new volume is: " << new_volume << std::endl;
   
}

void prey_stress_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    static Cell_Definition* pPreyDef = find_cell_definition( "Cell" );
    static Cell_Definition* pPredDef = find_cell_definition( "Predator" );
    
    // Increase stress levels if predator is nearby
    for( int n=0; n < pCell->state.neighbors.size(); n++ )
    {
        Cell* nbr = pCell->state.neighbors[n];
        if( nbr->type == pPredDef->type ){
            pCell -> custom_data["stress"] +=2;
//            std::cout << "Cell ID ("<< pCell->ID <<" ) detects predator. " << std::endl;
        }
    }
    
    // Increase stress levels if there isn't food around
    int food_index = microenvironment.find_density_index("food");
    double food_nearby = pCell->nearest_density_vector()[food_index];
    double threshold = 1.0;
    if (food_nearby < threshold ){
        pCell -> custom_data["stress"] +=0.1;
//        std::cout << "Cell ID ("<< pCell->ID <<" ) detects low food. " << std::endl;

    }
    else {
        pCell -> custom_data["stress"] -=0.5;
//        std::cout << "Cell ID ("<< pCell->ID <<" ) detects adequate food. " << std::endl;
    }
    
    
    // Update stress based on decay rate
    static double stress_decay_rate = pCell -> custom_data["stress_decay_rate"];
    pCell->custom_data["stress"] /= (1.0 + dt*stress_decay_rate);
    
    // Ensure stress does not go negative
    if (pCell->custom_data["stress"] < 0)
    {
        pCell->custom_data["stress"] = 0;
    }
//    std::cout << "Cell ID ("<< pCell->ID <<" ) new stress level is: " << pCell->custom_data["stress"] << std::endl;
    
}

Cell* perform_division(Cell* pCell){
    Cell* child = pCell -> divide(); // Perform division, creating a child cell
//    std::cout << "Division Occured! Parent Cell ID: " << pCell->ID << ". Child Cell ID: "<< child->ID <<std::endl;
    return child;
}

//
//void recursive_division_and_attach(PhysiCell::Cell* pCell, double initial_volume, int& division_count, std::vector<PhysiCell::Cell*>& lineage_cells) {
////    std::cout << "Evaluating Cell ID: " << pCell->ID << " for further division." << std::endl;
//
//    double current_volume = pCell->phenotype.volume.total;
//
//    // Check if the cell can divide
//    if (current_volume >= 2 * initial_volume) {
//        PhysiCell::Cell* child = perform_division(pCell);
////        std::cout << "Division Occurred! Parent Cell ID: " << pCell->ID << ". Child Cell ID: " << child->ID << std::endl;
//        division_count++;
//        
//        // Add the child to the lineage cells
//        lineage_cells.push_back(child);
//
//        // Attach child to all other cells in lineage
//        for(auto other_cell : lineage_cells) {
//            if(other_cell != child) {
//                child->attach_cell(other_cell);
////                std::cout << "Attaching Cell ID " << child->ID << " to Cell ID " << other_cell->ID << std::endl;
////                std::cout << "Cell (ID: " << child->ID << " ) has number of attachments =  " << child -> state.number_of_attached_cells() << std::endl;
//                
//            }
//        }
//        // Attach the parent to the child if not already done
//        pCell->attach_cell(child);
////        std::cout << "Attaching Cell ID " << pCell->ID << " to Cell ID " << child->ID << std::endl;
////        std::cout << "Cell (ID: " << pCell->ID << " ) has number of attachments =  " << pCell -> state.number_of_attached_cells() << std::endl;
//        
//        
//        // Recursive calls
//        recursive_division_and_attach(pCell, initial_volume, division_count, lineage_cells);
//        recursive_division_and_attach(child, initial_volume, division_count, lineage_cells);
//    }
//}



void iterative_division_and_attach(Cell* pCell, double initial_volume) {
    //    std::cout << "Evaluating Cell ID: " << pCell->ID << " for further division." << std::endl;
    std::vector<Cell*> current_cells;
    current_cells.push_back(pCell);
    
    int target_divisions = (int)floor(log2(pCell->phenotype.volume.total / initial_volume));
    
    
    for (int i = 0; i < target_divisions; ++i) {
        std::vector<Cell*> new_cells;
        
        // Perform division for each cell and update relationships
        for (auto cell : current_cells) {
            if (cell->phenotype.volume.total >= 2 * initial_volume) {
                Cell* child = perform_division(cell);
                new_cells.push_back(child);
                cell->attach_cell(child); // Attach parent to child
                child->attach_cell(cell); // Attach child to parent
                
                // If you are attached to others, label 'inedible'
                if (cell->state.number_of_attached_cells() > 0){
                    cell -> custom_data["edible"] = false;
                }
                else {cell -> custom_data["edible"] = true;}
                
                
                // Label child as inedible if attached
                if (child->state.number_of_attached_cells() > 0){
                    child -> custom_data["edible"] = false;
                }
                else {child -> custom_data["edible"] = true;}
                
                
                
                // Also, attach new child to all previous lineage cells
                for (auto lineage_cell : current_cells) {
                    if (lineage_cell != cell) { // No need to attach a cell to itself
                        child->attach_cell(lineage_cell);
                        lineage_cell->attach_cell(child);
                    }
                }
            }
            // If no division is needed, still carry over the cell to the next cycle
            new_cells.push_back(cell); // Update the list of current cells for next division cycle
        }
        current_cells = new_cells; // Update the list of current cells
    }
}



void mutation_function( Cell* pCell, Phenotype& phenotype, double dt ){
    
//    std::cout << "Mutation Function called for Cell ID: "<< pCell->ID << std::endl;
    
    // Mutate Transition rate:
    double multiplier2= 0.95 + 0.1*uniform_random();
    
    // Mutate life cycle gene
    pCell->custom_data["life_cycle_p"] *= multiplier2;
    // Clamp life_cycle_p to not exceed 1
    pCell->custom_data["life_cycle_p"] = std::min(1.0,pCell->custom_data["life_cycle_p"] );
    
    // Retrieve custom data from the specific cell
    double life_history_t = pCell->custom_data["life_history_t"];
    double life_cycle_p = pCell->custom_data["life_cycle_p"];
    
    
    // Set transition rates for phase 0 to 1 and phase 2 to 3
    double phase_0_duration = (1 - life_cycle_p) * life_history_t;
    double phase_2_duration = life_cycle_p * life_history_t;
    
    // Modify the cell cycle of the specific cell
    pCell->phenotype.cycle.data.transition_rate(0, 1) = 1.0 / phase_0_duration;
    pCell->phenotype.cycle.data.transition_rate(2, 3) = 1.0 / phase_2_duration;
}


// ~~~~ Cycle Phase Transitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool G0_arrest_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    // If the cell volume has not reached at least twice initial_volume, do not proceed to next phase.
    double current_volume = pCell ->phenotype.volume.total;
    double initial_volume = pCell -> custom_data["initial_volume"];
 
    if (current_volume < 2* initial_volume ){
        return true;
    }
    return false;
}
bool S_arrest_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    // access stress level
    double stress_level = pCell -> custom_data["stress"];
    
    // Get the threshold for arrest function
    double stress_threshold = pCell -> custom_data["cycle_arrest_stress_threshold"];
    if (stress_level > stress_threshold ){
        return true;
    }
    return false;
    
    
}


void phase_0_entry_function(Cell* pCell, Phenotype& phenotype, double dt)
{return;}

void phase_0_exit_function(Cell* pCell, Phenotype& phenotype, double dt)
{//   std::cout << "Phase 0 Exit Function called for Cell ID: "<< pCell->ID << std::endl;
    
    // Disable motility
    pCell -> phenotype.motility.is_motile = false;
    pCell -> phenotype.motility.migration_speed = 0.0;
}

void phase_1_entry_function(PhysiCell::Cell* pCell, PhysiCell::Phenotype& phenotype, double dt)
{
//    std::cout << "Phase 1 Entry Function called for Cell ID: " << pCell->ID << std::endl;

    // Disable motility
    pCell->phenotype.motility.is_motile = false;
    pCell->phenotype.motility.migration_speed = 0.0;
    
    // In preperation for division, mutate cell (so that all children have same genes)
    mutation_function( pCell, phenotype, dt ); // Mutate transition rates
    
    double initial_volume = pCell->custom_data["initial_volume"];

    // Start the iterative division process
    iterative_division_and_attach(pCell, initial_volume);
    
    
//    double volume_before_division = pCell ->phenotype.volume.total;
//    double initial_volume = pCell->custom_data["initial_volume"];
//    int division_count = 0;
//    std::vector<PhysiCell::Cell*> lineage_cells;  // This will hold all cells in the current lineage
//    lineage_cells.push_back(pCell);  // Start with the original cell
//
//    // Start the recursive division and attaching process
//    recursive_division_and_attach(pCell, initial_volume, division_count, lineage_cells);
//
////    std::cout << "For Cell ID: "<< pCell->ID << " , the volume before division was = " << volume_before_division << " , and the number of divisions this step = " << division_count << std::endl;
//    
//    // Clear the lineage cells vector after all divisions are complete
//    lineage_cells.clear();
////    std::cout << "Cleared lineage cells, total elements now: " << lineage_cells.size() << std::endl;

    return;
}


void phase_1_exit_function( Cell* pCell, Phenotype& phenotype, double dt )
{
//    std::cout << "Phase 1 Exit Function called for Cell ID: " << pCell->ID << std::endl;
    
    // Make sure all cells are their parent's initial size after division
    double initial_volume = pCell->custom_data["initial_volume"];
    pCell -> set_total_volume( initial_volume);
    
    // Mobilize cell:
    pCell -> phenotype.motility.is_motile = true;
    pCell -> phenotype.motility.migration_speed = 1.0;
return; }


void phase_2_entry_function(Cell* pCell, Phenotype& phenotype, double dt)
{ //std::cout << "Phase 2 Entry Function called for Cell ID: " << pCell->ID << std::endl;
    
    // Mobilize cell:
    pCell -> phenotype.motility.is_motile = true;
    pCell -> phenotype.motility.migration_speed = 1.0;
}


void phase_2_exit_function(Cell* pCell, Phenotype& phenotype, double dt)
{ //    std::cout << "Phase 2 Exit Function called for Cell ID: " << pCell->ID << std::endl;
    return;}

void phase_3_entry_function(Cell* pCell, Phenotype& phenotype, double dt)
{ //    std::cout << "Phase 3 Entry Function called for Cell ID: " << pCell->ID << std::endl;
    return;}

void phase_3_exit_function(Cell* pCell, Phenotype& phenotype, double dt)
{ //    std::cout << "Phase 3 Exit Function called for Cell ID: " << pCell->ID << std::endl;
  // Remove attachments
    pCell -> remove_all_attached_cells();
}



// ****************************************************************************
// ****************************************************************************
// PREDATOR FUNCTONS
// ****************************************************************************
// ****************************************************************************

void pred_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ){
    
    
    
    // Hunt and eat prey cells (and increase energy)
    pred_hunt_function(pCell, phenotype, dt);
    
    // Update energy based on metabolism
    static double metabolic_rate = pCell -> custom_data["metabolic_rate"];
    pCell->custom_data["energy"] /= (1.0 + dt*metabolic_rate);
    
    
    // Death based on starvation
    if( pCell->custom_data["energy"] < 1.0 )
    {
        pCell->lyse_cell();
        
        
    }
}

void pred_hunt_function( Cell* pCell, Phenotype& phenotype, double dt ){
    static Cell_Definition* pPreyDef = find_cell_definition( "Cell" );
    static Cell_Definition* pPredDef = find_cell_definition( "Predator" );
    
    bool is_hungry = (pCell->custom_data["energy"] < 40);
    
    
//    std::cout << "Predator (ID: " << pCell->ID  << ") number of neighbors = " <<   pCell->state.neighbors.size() << std::endl;

    
    for( int n=0; n < pCell->state.neighbors.size(); n++ )
    {
       
        
        Cell* nbr = pCell->state.neighbors[n];

        if( nbr->type == pPreyDef->type && nbr->custom_data["edible"] == true && is_hungry ){
            pCell -> ingest_cell(nbr);
            
            
            pCell -> custom_data["energy"] +=10;
            
//            std::cout << "Predator has ingested prey: " << nbr->ID << ". who had number of attachments =  " << nbr -> state.number_of_attached_cells() << std::endl;
            
            is_hungry = (pCell->custom_data["energy"] < 40);  // Re-evaluate hunger in case the condition changes due to energy increase


        }
        
    }
}

bool pred_repo_arrest_function( Cell* pCell, Phenotype& phenotype, double dt ){
    // If the cell volume has not reached at least twice initial_volume, do not proceed to next phase.
    double energy = pCell ->custom_data["energy"];
 
    if (energy < 20 ){
        return true;
    }
    return false;
    
}


