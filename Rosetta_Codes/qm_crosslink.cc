// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smturzo/qm_crosslinker_for_cyclic_peptides.cc
/// @brief An Application to define crosslinkers in cyclic peptide with Quantum Chemistry package (GAMESS)
/// @author smturzo (turzo.1@osu.edu)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/constraint_generator/util.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/simple_moves/SetTorsion.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
// core headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh> 
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/vector1.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/file/file_sys_util.hh>

// basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/EulerAngles.hh>

// c++ external header
#include <Eigen/SVD>
#include <Eigen/Dense>

//Headers for the Rosetta thread manager
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>

static basic::Tracer TR("apps.pilot.qm_crosslink");

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring::constraints;
using namespace basic::thread_manager;


//Options (ugh -- global variables):
OPT_KEY (String, pose_setup_file)
OPT_KEY (String, dofs_file)
OPT_KEY (String, output_weights_file)
OPT_KEY (Integer, pert_sample_number)
OPT_KEY (Integer, num_threads_qm)
OPT_KEY (Boolean, dump)
OPT_KEY (Boolean, debug_matrices)
OPT_KEY (Real, discard_sample_by_fa_repulsion_cutoff)
OPT_KEY (Real, delta_qm_energy_cutoff)

/// @brief This is an enum class. To ensure we are comparing integers instead of string.
/// It is efficient to do it this way.
/// MDIHEDRAL stands for Manual-Dihedral. This means user can set the dihedral angle for perturbation of any 4 allowed atoms. 
enum class PerturbedDofType { MAINCHAIN=1, SIDECHAIN, JUMP, BONDANGLE, MDIHEDRAL, BONDLENGTH, number_of_options=BONDLENGTH };

// This class is main Perturbation class. The class needs to be called with string and pose. It does a few things.
// It makes sure the dof file is correctly set up. Or else it will throw errors that will guide users to the right structure of the dof file
// This class also contains the perturb_pose member function that defines how each perturbation is supposed to happen. 
// This class may have room for adding other types of perturbation to the dof file.
// As well as defining how the new perturbations defined in dof file should behave in the member function (perturb_pose).
class Perturbation {

public:

	Perturbation() = delete; // Deleting default constructor.
	// Line Parsing constructor
	// To ensure strings get parsed once and when perturbation is applied only integers compared.
	// This is good coding practice. Improves the efficiency of the code by 50x or so.

    // Line parsing constructor.
    // This ensures that the string gets parsed once and once only,
    // and when perturbations are applied, you only compare integers.
    // All these accompanied by a number of checks! So needs to be double checked if they make sense.
    Perturbation( std::string const & line, core::pose::Pose const & my_pose ) {
      std::string const errmsg( "Error in Perturbation class constructor: " );
      std::istringstream iss( line );
      utility::vector1< std::string > linesplit = utility::split_whitespace(iss.str());
      std::string lineheader;
      core::Size const total_number_of_residues = my_pose.size();
      iss >> lineheader;
      if( lineheader == "MAINCHAIN" ) {
            if (linesplit.size() != 4) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For perturbation along main chain the DOF line should contain: Perturbation_Type(string, this should be: MAINCHAIN) Residue_Number(int) Dihedral_index(int) Angle_Perturbation(float).");
            dof_type_ = PerturbedDofType::MAINCHAIN;
            iss >> res_number_ >> dof_number_ >> pert_magnitude_;
            // Let's do some basic checks for the DOF File entered.
            check_residue_num_in_dof_file( errmsg, res_number_, total_number_of_residues, lineheader );
            core::Size total_torsion_number = my_pose.residue(res_number_).mainchain_torsions().size(); // Needs to be here because we want to check if the res numbers are valid first, before getting here.
            check_dof_num_in_dof_file( errmsg, res_number_, dof_number_, total_torsion_number, lineheader );

        } else if ( lineheader == "SIDECHAIN") {
            if (linesplit.size() != 4) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For perturbation along side chain the DOF line should contain: Perturbation_Type(string, this should be: SIDECHAIN) Residue_Number(int) Dihedral_index(int) Angle_Perturbation(float)."); 
            dof_type_ = PerturbedDofType::SIDECHAIN;
            iss >> res_number_ >> dof_number_ >> pert_magnitude_;
            // Let's do some basic checks for the DOF File entered.
            check_residue_num_in_dof_file( errmsg, res_number_, total_number_of_residues, lineheader );
            core::Size total_chi_torsion_number = my_pose.residue(res_number_).nchi(); // Needs to be here because we want to check if the res numbers are valid first, before getting here.
            check_dof_num_in_dof_file( errmsg, res_number_, dof_number_, total_chi_torsion_number, lineheader );

        } else if ( lineheader == "JUMP") {
            if (linesplit.size() != 4) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For perturbation along side chain the DOF line should contain: Perturbation_Type(string, this should be: JUMP) Jump#(int) Jump_Perturbation_Tranlation_Magnitude(float) Jump_Perturbation_Rotation_Magnitude(float)."); 
        	dof_type_ = PerturbedDofType::JUMP;
            iss >> jump_number_ >> jump_offset_tranlation_ >> jump_offset_rotation_;
            // Let's do some basic checks for the DOF File entered.
            check_jump_dof_in_dof_file( errmsg, jump_number_);

        } else if ( lineheader == "BONDANGLE") {
        	if (linesplit.size() != 8) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For bond angle perturbation, the DOF line should contain: Perturbation_Type(string, this should be: BONDANGLE) Residue_Number(int) Atom1_Name(string, atom type) Residue_Number(int) Atom2_Name(string, atom type) Residue_Number(int) Atom3_Name(string, atom type) Bond_Angle_Perturbation_Magnitude(float)");
        	dof_type_ = PerturbedDofType::BONDANGLE;
        	iss >> res_number_ >> atom_name_ >> res2_number_ >> atom2_name_ >> res3_number_ >> atom3_name_ >> pert_magnitude_;
        	// Let's do some basic checks for the DOF File entered.
            check_residue_num_in_dof_file( errmsg, res_number_,  total_number_of_residues, lineheader );
            check_residue_num_in_dof_file( errmsg, res2_number_, total_number_of_residues, lineheader );
            check_residue_num_in_dof_file( errmsg, res3_number_, total_number_of_residues, lineheader );                        
            // if atom names does not exist in residue utility exist out
            atom_number_  = check_atom_name_in_dof_file_set_atom_num( errmsg, res_number_,  atom_name_,  lineheader, my_pose );
            atom2_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res2_number_, atom2_name_, lineheader, my_pose );
            atom3_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res3_number_, atom3_name_, lineheader, my_pose );

        } else if (lineheader =="MDIHEDRAL") {
            if (linesplit.size() != 10) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For maunual dihedral perturbation, the DOF line should contain: Perturbation_Type(string, this should be: MDIHDERAL) Residue_Number(int) Atom1_Name(string, atom type) Residue_Number(int) Atom2_Name(string, atom type) Residue_Number(int) Atom3_Name(string, atom type) Residue_Number(int) Atom4_Name(string, atom type) Bond_Angle_Perturbation_Magnitude(float)");
            dof_type_ = PerturbedDofType::MDIHEDRAL;
            iss >> res_number_ >> atom_name_ >> res2_number_ >> atom2_name_ >> res3_number_ >> atom3_name_ >> res4_number_ >> atom4_name_ >> pert_magnitude_;
            // Let's do some basic checks for the DOF File entered.
            check_residue_num_in_dof_file( errmsg, res_number_,  total_number_of_residues, lineheader );
            check_residue_num_in_dof_file( errmsg, res2_number_, total_number_of_residues, lineheader );
            check_residue_num_in_dof_file( errmsg, res3_number_, total_number_of_residues, lineheader );                        
            check_residue_num_in_dof_file( errmsg, res4_number_, total_number_of_residues, lineheader );                        
            // if atom names does not exist in residue utility exist out
            atom_number_  = check_atom_name_in_dof_file_set_atom_num( errmsg, res_number_,  atom_name_,  lineheader, my_pose );
            atom2_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res2_number_, atom2_name_, lineheader, my_pose );
            atom3_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res3_number_, atom3_name_, lineheader, my_pose );
            atom4_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res4_number_, atom4_name_, lineheader, my_pose );

        } else if ( lineheader == "BONDLENGTH") {
        	if (linesplit.size() != 6) utility_exit_with_message(errmsg+"Dof file setup is wrong. User error. For bond length perturbation, the DOF line should contain: Perturbation_Type(string, BONDLENGTH) Residue_Number(int) Atom1_Name(string, atom type) Residue_Number(int) Atom2_Name(string, atom type) Bond_Length_Perturbation_Magnitude(float, this cannot be negative).");
        	dof_type_ = PerturbedDofType::BONDLENGTH;
        	iss >> res_number_ >> atom_name_ >> res2_number_ >> atom2_name_ >> pert_magnitude_;
        	// Let's do some basic checks for the DOF File entered.
            check_residue_num_in_dof_file( errmsg, res_number_,  total_number_of_residues, lineheader );
            check_residue_num_in_dof_file( errmsg, res2_number_, total_number_of_residues, lineheader );
            atom_number_  = check_atom_name_in_dof_file_set_atom_num( errmsg, res_number_,  atom_name_,  lineheader, my_pose );
            atom2_number_ = check_atom_name_in_dof_file_set_atom_num( errmsg, res2_number_, atom2_name_, lineheader, my_pose );

        } else utility_exit_with_message( errmsg + "Did not recognize perturbation type " + linesplit[1] + "!" );
    }	

public:
    // Fun stuff here: perturb pose 
    // Function that actually does the perturbation.  This is called
    // hundreds or thousands or millions of times, so it shouldn't do
    // any string parsing!
    void perturb_pose( core::pose::Pose & my_pose ) const {
      core::Real current_angle(0);
      core::Real new_angle(0);
			switch( dof_type_ ) {  
				case PerturbedDofType::MAINCHAIN :
					current_angle = my_pose.Pose::torsion( core::id::TorsionID( res_number_, core::id::BB, dof_number_ ) );
					new_angle = current_angle+ (numeric::random::gaussian()*pert_magnitude_);
					my_pose.set_torsion(core::id::TorsionID( res_number_, core::id::BB, dof_number_ ),new_angle);
					break;

				case PerturbedDofType::SIDECHAIN :
					current_angle = my_pose.Pose::torsion( core::id::TorsionID( res_number_, core::id::CHI, dof_number_ ) );
					new_angle = current_angle + (numeric::random::gaussian()*pert_magnitude_);
					my_pose.set_torsion(core::id::TorsionID( res_number_, core::id::CHI, dof_number_ ),new_angle); // double check degrees
					break;

				case PerturbedDofType::JUMP : {
					core::kinematics::Jump jump_pose( my_pose.jump(jump_number_) );
					numeric::xyzVector < core::Real > transvect = jump_pose.get_translation();
					transvect[1] += ( numeric::random::gaussian()* jump_offset_tranlation_ );
					transvect[2] += ( numeric::random::gaussian()* jump_offset_tranlation_ );
					transvect[3] += ( numeric::random::gaussian()* jump_offset_tranlation_ );
					core::Real const jump_offset_rotation_radians( numeric::conversions::radians(jump_offset_rotation_) ); 
					numeric::EulerAngles < core::Real > euler_angles( jump_pose.get_rotation() );
					euler_angles.phi( numeric::principal_angle_radians( euler_angles.phi() + numeric::random::gaussian() * jump_offset_rotation_radians ) ); 
					euler_angles.psi( numeric::principal_angle_radians( euler_angles.psi() + numeric::random::gaussian() * jump_offset_rotation_radians ) ); 
					euler_angles.theta( std::min( numeric::constants::d::pi, std::max( 0.0, euler_angles.theta() + numeric::random::gaussian() * jump_offset_rotation_radians ) ) ); 
					jump_pose.set_translation( transvect ); // set the new translation to the copied jump
					jump_pose.set_rotation( euler_angles.to_rotation_matrix() ); // set the new rotation to the copied jump
					my_pose.set_jump( jump_number_, jump_pose ); // Perturb the pose Jump
					break;
				}
				
				case PerturbedDofType::BONDANGLE : {
						core::Real const pert_magnitude_radians( numeric::conversions::radians(pert_magnitude_) ); 
						core::id::AtomID atom_id1( core::id::AtomID( atom_number_,  res_number_  ) );
						core::id::AtomID atom_id2( core::id::AtomID( atom2_number_, res2_number_ ) );
						core::id::AtomID atom_id3( core::id::AtomID( atom3_number_, res3_number_ ) );
						core::Real current_angle = my_pose.conformation().bond_angle( atom_id1, atom_id2, atom_id3);
						core::Real new_angle = current_angle + ( numeric::random::gaussian()*pert_magnitude_radians ) ; // new angle in radians
						my_pose.conformation().set_bond_angle( atom_id1, atom_id2, atom_id3, new_angle ); // Set new bond angle (radians).
						break;
				}

                case PerturbedDofType::MDIHEDRAL : {
                        core::Real const pert_magnitude_radians( numeric::conversions::radians( pert_magnitude_ ) ); 
                        core::id::AtomID atom_id1( core::id::AtomID( atom_number_,  res_number_  ) );
                        core::id::AtomID atom_id2( core::id::AtomID( atom2_number_, res2_number_ ) );
                        core::id::AtomID atom_id3( core::id::AtomID( atom3_number_, res3_number_ ) );
                        core::id::AtomID atom_id4( core::id::AtomID( atom4_number_, res4_number_ ) );
                        current_angle = my_pose.conformation().torsion_angle( atom_id1, atom_id2, atom_id3, atom_id4 ); // find torsion angle for a given set of atoms manually
                        new_angle = current_angle + ( numeric::random::gaussian()*pert_magnitude_radians ) ; // new angle in radians
                        my_pose.conformation().set_torsion_angle( atom_id1, atom_id2, atom_id3, atom_id4, new_angle ); // Set new dihedral angle (radians).
                        break;
                }

				case PerturbedDofType::BONDLENGTH : {
                        core::id::AtomID atom_id1( core::id::AtomID( atom_number_,  res_number_  ) );
                        core::id::AtomID atom_id2( core::id::AtomID( atom2_number_, res2_number_ ) );
                        
                        core::Real current_bond_length( my_pose.conformation().bond_length( atom_id1, atom_id2 ) );
						//core::Real current_bond_length(my_pose.residue( res_number_ ).xyz( atom_number_ ).distance( my_pose.residue( res2_number_ ).xyz( atom2_number_ ))); // get current bond length
						core::Real new_bond_length = std::abs( current_bond_length + ( numeric::random::gaussian()*pert_magnitude_ ) ); // get new bond length. std::abs avoids the negative bond length.
						my_pose.conformation().set_bond_length( atom_id1, atom_id2, new_bond_length ); // Set new bond length. TODO: This doesnt work for atoms that are not bonded. Need to add checks to utility exit with message if users such atoms are added. Add this check to perturbation list class.
						break;
				}    
		}
		my_pose.update_residue_neighbors();
  } 

private:

    void check_residue_num_in_dof_file( std::string const &error_message_to_throw, core::Size const residue_number, core::Size const pose_tot_residues, std::string const &dof_lineheader ) {
        if ( residue_number >  pose_tot_residues ) utility_exit_with_message( error_message_to_throw+ "Residue number "+ std::to_string( residue_number ) + " in the line " + dof_lineheader + " in dof file cannot be greater than the size of the sequence length" );
        if ( res_number_<= 0 ) utility_exit_with_message( error_message_to_throw + "Residue number " + std::to_string( residue_number ) + " in the line " + dof_lineheader + " in dof file cannot be equal or less than zero." );
        if ( res_number_< 1 ) utility_exit_with_message( error_message_to_throw + "Residue number " + std::to_string( residue_number ) + " in the line "+ dof_lineheader + " in dof file cannot be less than 1." );            

    }
    core::Size check_atom_name_in_dof_file_set_atom_num( std::string const &error_message_to_throw, core::Size const residue_number, std::string const &atom_name, std::string const &dof_lineheader, core::pose::Pose const & my_pose ) {
        if (my_pose.residue( residue_number ).has( atom_name ) ==1) {
            core::Size atom_number(my_pose.residue_type( residue_number ).atom_index( atom_name ) );
            return atom_number;
        } else utility_exit_with_message( error_message_to_throw +"Atom name "+atom_name+" in the line "+dof_lineheader+" in the dof file could not be found in pose. In a pdb file these are usually defined in the column right next to the atom numbering.");
    }

    void check_dof_num_in_dof_file( std::string const &error_message_to_throw, core::Size const residue_number, core::Size const dof_number, core::Size const total_torsion_to_perturb, std::string const &dof_lineheader ) {
        if ( dof_number <=0 ) utility_exit_with_message( error_message_to_throw + "Torsion index cannot be equal or less than zero in dof file. User error. Torsion index  residue " +std::to_string(residue_number)+" in the line "+dof_lineheader+" needs be an integer between 1 <= torsion_index <= "+ std::to_string( total_torsion_to_perturb ) );
        if ( dof_number <1 ) utility_exit_with_message( error_message_to_throw + "Torsion index cannot be less than 1 in dof file. User error. Torsion index of residue " + std::to_string(residue_number) + " in the line " + dof_lineheader + " needs be an integer between 1 <= torsion_index <= " + std::to_string( total_torsion_to_perturb ) );
        if ( dof_number > total_torsion_to_perturb ) utility_exit_with_message( error_message_to_throw + "Torsion index cannot be greater than the total number of torsions in residue number " + std::to_string(residue_number) + " in the line " + dof_lineheader + ". User error. Total number of torsion index is " + std::to_string(total_torsion_to_perturb ) );
    }

    void check_jump_dof_in_dof_file( std::string const &error_message_to_throw, core::Size const jump_number ) {
        if ( jump_number < 1 ) utility_exit_with_message( error_message_to_throw+ "Jump perturbation is indicated in dof file, but jump number is set to  "+ std::to_string(jump_number)+ "in dof file. User error. Jump number has to be >= 1." );
        if ( jump_number <= 0 ) utility_exit_with_message( error_message_to_throw+ "Jump perturbation is indicated in dof file, but jump number "+ std::to_string(jump_number)+ "in dof file is negative. User error. Jump number has to be >= 1." ); 
    }

private:

	PerturbedDofType dof_type_ = PerturbedDofType::MAINCHAIN;
	core::Size dof_number_ = 0;
	core::Size res_number_ = 0;
	core::Size res2_number_ = 0;
    core::Size res3_number_ = 0;
	core::Size res4_number_ = 0;

	std::string atom_name_  = "";
	core::Size atom_number_  = 0;
	std::string atom2_name_ = "";
	core::Size atom2_number_ = 0;
	std::string atom3_name_ = "";
	core::Size atom3_number_ = 0;
    std::string atom4_name_ = "";
    core::Size atom4_number_ = 0;

	core::Real pert_magnitude_ = 0.0;
	core::Size jump_number_ = 0;
	core::Real jump_offset_tranlation_ = 0.0;
	core::Real jump_offset_rotation_ = 0.0;
};

// Is it Rosetta convention to put typedef at the top of the file ?
typedef utility::pointer::shared_ptr< Perturbation > PerturbationOP;
typedef utility::pointer::shared_ptr< Perturbation const > PerturbationCOP;

// This class needs to be called with a vector of strings and pose. 
// Within this class the member function will then perturb the pose one by one.
class PerturbationList
{
public:
    // Delete the default constructor.  Forces devs to use the
    // constructor that parses lines.
    PerturbationList() = delete;

    // Constructor that parses the contents of a perturbation file.
    // This calls the line-parsing constructor of the Perturbation class when it creates
    // each perturbation.  This ensures that the Perturbation parses the strings once and only
    // once.  The pose is used to ensure that the perturbed DoF actually exists, and to throw an
    // informative error if it does not.
    PerturbationList( utility::vector1< std::string > const & lines, core::pose::Pose const & my_pose ) {
        for( std::string const & line : lines ) {
            perturbations_.push_back( utility::pointer::make_shared<Perturbation>(line, my_pose) );
        }
    }

public:

    // Function that applies all the perturbations.  This is called many times.
    // It calls Perturbation::perturb_pose(), which does no string-parsing.
    void perturb_pose( core::pose::Pose & my_pose ) {
        for( PerturbationOP const & pert : perturbations_ ) {
            pert->perturb_pose( my_pose );
        }
    }

public:


private:
    utility::vector1< PerturbationOP > perturbations_;
};

// This class has to be initialized with the number of samples (int) and the constraint set
// This class has 3 member functions: do_sampling(...), lin_alg_least_square_fitting(), write_report_from_sampling(...)
// In this class, we could potentially put a future member function that would allow generation of a lot of data for Machine Learning
class PerturbationSampler {

public:
	PerturbationSampler() = delete;
	PerturbationSampler(core::Size const number_pert_samples, ConstraintSetCOP cst) :
	number_pert_samples_(number_pert_samples), 
    cst_(cst), 
    cst_list_( cst_->get_all_constraints() ), 
    cst_size_(cst_list_.size()), 
    sampled_matrix_(number_pert_samples_,(cst_size_)), 
    qm_vector_(number_pert_samples_),
    constraint_weights_vector_( (cst_size_+1) )
	{

	}

public:
	void calculate_qm_score_in_thread(
        core::pose::Pose const & original_pose, 
        core::scoring::ScoreFunction const & original_sfxn, 
        PerturbationList & perturbation, 
        core::Size const this_job_index, 
        Eigen::VectorXd & scores_out,  
        core::Real fa_repulsion_score_cuttoff, 
        bool const do_dump_pdb, 
        Eigen::MatrixXd & sampled_matrix, 
        utility::vector0<bool> & samples_successful) const {
    		core::pose::PoseOP pose_copy( utility::pointer::make_shared< core::pose::Pose >() ); //Made an empty pose.
    		pose_copy->detached_copy( original_pose );
    		//Perturb pose_copy here.
    		perturbation.perturb_pose( *pose_copy );
    		//if (fa_rep_cutoff_bool) {
        	 // Created an empty scorefunction, then set fa_rep score to 1.
            core::scoring::ScoreFunctionOP fa_rep_sfxn(utility::pointer::make_shared< core::scoring::ScoreFunction >());
            fa_rep_sfxn->set_weight(core::scoring::score_type_from_name("fa_rep"), 1.0 );
            core::Real fa_rep_score( (*fa_rep_sfxn)(*pose_copy) );
            TR.Info << "The FA repulsion score of the sample is: " << fa_rep_score << std::endl;
            samples_successful[this_job_index-1] = fa_rep_score < fa_repulsion_score_cuttoff ;
            //}

            if (do_dump_pdb) pose_copy->dump_pdb("PERTURUBED_SAMPLE_"+std::to_string(this_job_index)+".pdb"); // ### ATTN WHY IS THIS ONLY OUTPUTTING THE LAST STRUCTURE?????   
    		// Score the pose with QM:
    		core::scoring::ScoreFunctionOP sfxn_copy( original_sfxn.clone() );
    		scores_out(this_job_index-1) = (*sfxn_copy)(*pose_copy);
            // Start a new constr_sfxn in each thread
            core::scoring::ScoreFunctionOP constr_sfxn(utility::pointer::make_shared< core::scoring::ScoreFunction >());
            constr_sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
            constr_sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
            constr_sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );        
    		// - Pass in the matrix of constraint scores.
            TR.Info << "Does it break here 1" << std::endl;
            for (core::Size cstr = 1; cstr <= cst_size_; cstr++) {
        		pose_copy->add_constraint( cst_list_[cstr] );
                core::Real constr_score( (*constr_sfxn)(*pose_copy) );
                sampled_matrix( (this_job_index-1), (cstr-1) ) = (constr_score);
            }
            TR.Info << "Does it break here 2" << std::endl;
    //sampled_matrix_( (this_job_index-1), (cst_size_) ) = 1.0; Don't need this here
	}

    /// @brief 
    /// @param my_pose 
    /// @param pert_list_object 
    /// @param qm_sfxn 
    /// @param fa_repulsion_score_cuttoff 
    /// @param do_dump_pdb 
    /// @param write_all_matrices 
    /// @param user_specified_delta_qm_energy 
    /// @param num_threads_to_use 
    void do_sampling(
        core::pose::Pose const & my_pose, 
        PerturbationList & pert_list_object, 
        core::scoring::ScoreFunctionOP qm_sfxn, 
        core::Real const fa_repulsion_score_cuttoff, 
        bool const do_dump_pdb,  
        bool const write_all_matrices, 
        core::Real const user_specified_delta_qm_energy, 
        core::Size const num_threads_to_use) {
        // Assignment info since it is application. Is this correct? What does it mean? What is it doing?
        RosettaThreadAssignmentInfo assignment_info( RosettaThreadRequestOriginatingLevel::APPLICATIONS_OR_APPLICATION_PROTOCOLS );
        utility::vector1< RosettaThreadFunction > work_vector; // Initiated a work vector for Rosetta Threading
        utility::vector0< bool > samples_successful( number_pert_samples_, false );
        if (has_sampling_been_done_) utility_exit_with_message("Sampling has been done. Exiting."); 
        //Starting for_loop for the total number of samples given by user.
        TR.Info << "Looping through Initial NSAMP: " << std::endl;
        for ( core::Size nsamp=1; nsamp<=number_pert_samples_; nsamp++ ) { 
            //The ampersand means "memory address of the function".
            //Non-static class member functions need to know which instance of the class they're being called from.
            work_vector.emplace_back( 
                std::bind(
                    &PerturbationSampler::calculate_qm_score_in_thread, 
                    this, 
                    std::cref( my_pose ), 
                    std::cref( *qm_sfxn ), 
                    std::ref ( pert_list_object ), 
                    nsamp, 
                    std::ref ( qm_vector_ ),  
                    fa_repulsion_score_cuttoff, 
                    do_dump_pdb, 
                    std::ref ( sampled_matrix_ ), 
                    std::ref ( samples_successful ) 
                ) 
            );
        }
        basic::thread_manager::RosettaThreadManager::get_instance()->do_work_vector_in_threads( work_vector, num_threads_to_use, assignment_info ); // Won't move on until all threads have finished their work.
        // At this point, all the threaded work is done, and the results are in the qm_vector_ vector (and will be in the sampled_matrix_ too.)
        // Before moving on, we should delete the rows that correspond to failed samples, based on the samples_successful vector.
        // TODO: Need to double check this. Vikram please help.
        // Heavily commented for clarification
        // We initiate a vector1 of integers called "indices_true"
        // In this vector we store all the indices of the vector0 samples_successful if the corresponding element of the index was True
        TR.Info << "Size of Sampled Matrix: Rows: " << sampled_matrix_.rows() << ", Columns: " << sampled_matrix_.cols() << std::endl;
        TR.Info << "Size of QM Vector: " << qm_vector_.size() << std::endl;
        TR.Info << "Size of Samples Successful: " << samples_successful.size() << std::endl; // This checksou okay!
        utility::vector1<core::Size> indices_true;
        for (core::Size count=1; count<=samples_successful.size(); count++ ) {
            TR.Info << "Sample Count " << count << " is " << samples_successful[count-1] << std::endl;
            if (samples_successful[count-1]) { // This should be fixed now, however I don't like the default value of 9999999999999999999999999999.9 for cutoff score. Is there a better way to do this??
                indices_true.push_back(count-1);
            } 
        }
        // Count the total number of samples with fa_rep cutoff
        core::Size count_successful_sample_true(indices_true.size());
        TR.Info << "SEE INDICES_TRUE VECTOR BELOW" << std::endl;
        TR.Info << "SEE INDICES_TRUE VECTOR BELOW" << std::endl;
        TR.Info << "SEE INDICES_TRUE VECTOR BELOW" << std::endl;
        for (core::Size idx_true=1; idx_true<= indices_true.size();idx_true++) {
            TR.Info << "IS Indices_True Going Negative? : " << indices_true[idx_true] << std::endl;
        }
        TR.Info << "SEE INDICES_TRUE VECTOR ABOVE" << std::endl; 
        TR.Info << "SEE INDICES_TRUE VECTOR ABOVE" << std::endl; 
        TR.Info << "SEE INDICES_TRUE VECTOR ABOVE" << std::endl;

        TR.Info << "Number of samples that passed fa_rep: " << count_successful_sample_true << std::endl;
        // Initiate Eigen vector and Matrix with the new number of rows (column stays the same) based on number of samples with fa_rep_cutoff (i.e. count_successful_sample_true variable)
        Eigen::VectorXd fa_rep_reduced_qm_vector(count_successful_sample_true); // row length is reduced because of fa_rep
        Eigen::MatrixXd fa_rep_reduced_samp_matrix(count_successful_sample_true,cst_size_); // row length is reduced because of fa_rep and column length one less (since we are not doing linear alg just yet)
        TR.Info << "Initiated Matrix fa_rep_reduced_samp_matrix size is: Rows: " << fa_rep_reduced_samp_matrix.rows() << ", Columns: "<< fa_rep_reduced_samp_matrix.cols() << std::endl;
        TR.Info << "Initiated Vector fa_rep_reduced_qm_vector size is: " << fa_rep_reduced_qm_vector.size() << std::endl;
        // Now we iterate over the total number of fa rep reduced samples
        // and to the new fa_rep_reduced_qm_vector we add the qm energy where index was true from the original qm vector (qm_vector_)
        for (core::Size reduced_samp = 1; reduced_samp <= count_successful_sample_true; reduced_samp++) {
            //if (qm_vector_(indices_true[reduced_samp]) < threshold_value) {}
            fa_rep_reduced_qm_vector(reduced_samp-1) = qm_vector_(indices_true[reduced_samp]);
            for (core::Size cstr = 1; cstr <= cst_size_; cstr++) {
                // similarly to the new fa_rep_reduced_samp_matrix we add the constr energy where index was true from the original qm vector (qm_vector_)
                fa_rep_reduced_samp_matrix((reduced_samp-1),(cstr-1)) = sampled_matrix_( (indices_true[reduced_samp]), (cstr-1) ); // TODO:: Vikram can you check if these indices make sense.    
            }
            //fa_rep_reduced_samp_matrix( (reduced_samp-1), (cst_size_)) = 1.0; // We should not be doing this here.
        }
        TR.Info << "fa_rep_reduced_samp_matrix is filled (?): " << fa_rep_reduced_samp_matrix.all() << std::endl;
        TR.Info << "fa_rep_reduced_qm_vector is filled (?): " << fa_rep_reduced_qm_vector.all() << std::endl;

        // Now find the lowest qm energy of the fa_rep_reduced_qm_vector
        // And set threshold value to lowest energy found plus tolerance defined by user.
        core::Real const lowest_enrg_found = fa_rep_reduced_qm_vector.minCoeff();
        core::Real const threshold_value = lowest_enrg_found + user_specified_delta_qm_energy;
    
        // We follow the same approach as before, i.e. We initiate a vector1 of integers called "qm_threshold_indices"
        // To this vector we add all indices where qm energy is less than a threshold value
        utility::vector1<core::Size> qm_threshold_indices;
        for (core::Size qm_threshold=1; qm_threshold<=count_successful_sample_true; qm_threshold++ ) {
            if (fa_rep_reduced_qm_vector(qm_threshold-1) < threshold_value)   {
                qm_threshold_indices.push_back(qm_threshold-1);
            } 
        }
        // Next we obtain a reduced number of samples based on qm energy threshold value
        core::Size count_sample_by_qmthreshold(qm_threshold_indices.size());

        // Using the new reduced number of samples we again create a new Eigen Vector and Matrix. (that is reduced number of rows by QM energy threshold)
        Eigen::VectorXd new_qm_vector(count_sample_by_qmthreshold); 
        Eigen::MatrixXd new_sample_matrix(count_sample_by_qmthreshold,cst_size_+1); // Note now we can finally make a new sample matrix with one extra column and fill that with 1 (i.e. the offset to do least square fit)
        TR.Info << "Initiated Matrix new_sample_matrix size is: Rows: " << new_sample_matrix.rows() << ", Columns: "<< new_sample_matrix.cols() << std::endl;
        TR.Info << "Initiated Vector new_qm_vector size is: " << new_qm_vector.size() << std::endl;
        // Now we fill new_qm_vector and new_sample_matrix with elements corresponding to indices in the qm_threshold_indices vector1. 
        for (core::Size qm_reduced_samp = 1; qm_reduced_samp <= count_sample_by_qmthreshold; qm_reduced_samp++) {
            new_qm_vector(qm_reduced_samp-1) = fa_rep_reduced_qm_vector(qm_threshold_indices[qm_reduced_samp]);
            for (core::Size cstr = 1; cstr <= cst_size_; cstr++) {
                new_sample_matrix((qm_reduced_samp-1),(cstr-1)) = fa_rep_reduced_samp_matrix( (qm_threshold_indices[qm_reduced_samp]), (cstr-1) );    
            }
            new_sample_matrix( (qm_reduced_samp-1), (cst_size_)) = 1.0; // It finally makes sense to do this here, since there are no row reductions.
        }
        TR.Info << "new_sample_matrix is filled (?): " << new_sample_matrix.all() << std::endl;
        TR.Info << "new_qm_vector is filled (?): " << new_qm_vector.all() << std::endl;

        sampled_matrix_ = new_sample_matrix;//Throw away points that are above threshold value (total rows is same as new qm vector, n columns)
        qm_vector_ = new_qm_vector;//New QM Vector has rows as sampled matrix , 1 column)

        if (write_all_matrices) {
            utility::io::ozstream outstream;
            outstream.open("Debug_Sample_Matrix_before_resize");
            outstream << sampled_matrix_ << std::endl;
            outstream.close();
            outstream.open("Debug_QM_Vector_before_resize");
            outstream << qm_vector_ << std::endl;
            outstream.close();
        }
        if (fa_repulsion_score_cuttoff > 0.0) TR.Info << "User opted for the option to reduce samples by fa_rep cutoff. Now the total number of samples is reduced to "+std::to_string(count_successful_sample_true)+" from original number of samples of "+ std::to_string(number_pert_samples_) << std::endl;
        TR.Info << "User opted for the option to refine QM energy by user defined delta QM energy"+std::to_string(user_specified_delta_qm_energy)+". Now the total number of samples is reduced to "+std::to_string(count_sample_by_qmthreshold)+" from original number of samples of "+ std::to_string(number_pert_samples_) << std::endl;
        //sampled_matrix_.conservativeResize(count_sample_by_qmthreshold,cst_size_+1); // Resized after threshold was applied. Dont need this.
        //qm_vector_.conservativeResize(counter_reduced_samples_by_threshold); // to get rid of bad sample for both the matrix and vector. Dont need this.
        has_sampling_been_done_ =true;
        if (write_all_matrices) {
            utility::io::ozstream outstream;
            outstream.open("Debug_Sample_Matrix_after_resize");
            outstream << sampled_matrix_ << std::endl;
            outstream.close();
            outstream.open("Debug_QM_Vector_after_resize");
            outstream << qm_vector_ << std::endl;
            outstream.close();
        }
    }
	// This is a separate function in case the devs want to try other function approximator (?)
	void lin_alg_least_square_fitting(bool const write_all_matrices) {
		if (!has_sampling_been_done_) utility_exit_with_message("Sampling has not been done. Exiting.");
		Eigen::BDCSVD<Eigen::MatrixXd> svd(sampled_matrix_,Eigen::ComputeFullU | Eigen::ComputeFullV);
        //constraint_weights_vector_.conservativeResize( qm_vector_.size() );
		constraint_weights_vector_ = svd.solve(qm_vector_);
		
		has_linear_algebra_been_done_=true;
        if (write_all_matrices) {
            utility::io::ozstream outstream;
            outstream.open("Debug_Constraint_weights_vector");
            outstream << constraint_weights_vector_ << std::endl;
            outstream.close();
        }
	}

	void write_report_from_sampling(std::string const & out_file_name) {
		if (!has_sampling_been_done_) utility_exit_with_message("Sampling has not been done. Exiting.");
		if (!has_linear_algebra_been_done_) utility_exit_with_message("Linear algebra has not been done. Will not write report.");
		utility::io::ozstream outstream;
		outstream.open( out_file_name );
		outstream << "CONSTRAINT_NUMBER: OPTIMIZED_WEIGHT" << std::endl;
        core::Size constr_vec_size( constraint_weights_vector_.size() );
		for( core::Size wts = 1; wts <= ( constr_vec_size - 1 ); wts++ ) {
            TR << "Writing Optimized weights for the " << wts << " of the cst file from constraing weight vector: " << (wts-1) << std::endl;
			outstream <<"CONSTRAINT_"<< wts << ": " << constraint_weights_vector_(wts-1) << std::endl;
		}
        TR << "Writing QM Constant after the optimized weights " << (constr_vec_size) << "  " << (constr_vec_size-1) << std::endl;
		outstream <<"QM_CONSTANT_"<< ( constr_vec_size ) << ": " << constraint_weights_vector_( constr_vec_size -1 ) << std::endl;
		outstream.close();
	}

private:
	core::Size const number_pert_samples_; 
	ConstraintSetCOP cst_;
	utility::vector1< ConstraintCOP > const cst_list_;
	core::Size cst_size_;
	Eigen::MatrixXd sampled_matrix_;
	Eigen::VectorXd qm_vector_;
	Eigen::VectorXd constraint_weights_vector_;
	bool has_sampling_been_done_ = false;
	bool has_linear_algebra_been_done_ = false;

};

//This function simply registers all options for user input
void register_options() {
    option.add_relevant( score::weights ); 
    option.add_relevant( constraints::cst_file );
    option.add_relevant( in::file::fold_tree ); // Should I take this out, as it is no longer relevant but optional ?

	NEW_OPT(pose_setup_file, "A file containing an XML tag setting up a PeptideStubMover to build the pose. Required input.", "" );
	NEW_OPT(dofs_file, "A file containing the degrees of freedom for sampling the pose. Required input.", "" ); 
	NEW_OPT(pert_sample_number, "Integer (n) number indicating the total number of perturbation to sample. For example for a given pose_setup_file, fold_tree  cst_file, qm weights file, and dof file it will randomly perturb n many times and get constraint scores and one qm score.",100);
	NEW_OPT( dump, "Dumps PDB for every time the sample is perturbed. WARNING: This is costly and is only to be used for debugging purposes.", false );
    NEW_OPT( debug_matrices, "Writes out the sampled matrix, qm_vector, and the constraint_weights_vector to their individual files. Only to be used for debugging",false);
	NEW_OPT(output_weights_file, "An output file name to store the optimized constraint weight. If none is provide this is automatically set to OPTIMIZED_CONSTRAINT_WEIGHTS.out. The last line in this file is the QM constant.","OPTIMIZED_CONSTRAINT_WEIGHTS.out");
    NEW_OPT(discard_sample_by_fa_repulsion_cutoff, "User may find it useful to sample structure with fa_rep_cuttoff score, where low fa_rep term corresponds to low steric clash. One reason to do this is to make sure that QM package is not crashing or moving atom in a completely different place. By default this value is set to 0.0.",0.0);
    NEW_OPT(delta_qm_energy_cutoff, "User may find it useful to sample structure with a QM delta energy cuttoff, where a large cutoff indicates sampling perturbation across a massive energy landscape. A small cuttoff indicates sampling close to an energy minima. By default this is set to 0.0.",0.0);
    NEW_OPT(num_threads_qm, "User may want to parallelize the QM sampling process over multiple threads. This options allows it. If none is passed the value that gets passed is 0. Meaning it is going to search for all available threads. Note this is multi-threading on a single node and not MPI.",0);
}

// This function builds pose from a parsed xml protocol
void build_initial_pose_from_xml_tag(std::string &xml_tag_file, std::string const &errmsg,core::pose::Pose & my_pose) {
	if (!utility::file::file_exists(xml_tag_file)) utility_exit_with_message(errmsg+"The file: "+xml_tag_file+" passed with the option -pose_setup_file does not exist in path.");
	utility::vector1< std::string > xml_tag_file_name_white_spacesplit_vector = utility::split_whitespace(xml_tag_file);
	if (xml_tag_file_name_white_spacesplit_vector.size()==0) utility_exit_with_message(errmsg+"The application only found whitespaces given with the option -pose_setup_file so it will quit.");
	std::string xml_tag_file_content;
	utility::io::izstream xml_tag_file_stream( xml_tag_file );
	utility::slurp( xml_tag_file_stream, xml_tag_file_content );
	utility::vector1< std::string > xml_tag_file_content_whitespace_split_vector = utility::split_whitespace(xml_tag_file_content);
	if (xml_tag_file_content_whitespace_split_vector.size()==0) utility_exit_with_message(errmsg+"Mesa your XML tag file is empty");
	basic::datacache::DataMap data; // this is an empty datamap
	utility::tag::TagOP mytag(utility::pointer::make_shared< utility::tag::Tag >());
	mytag->read(xml_tag_file_content);
	protocols::rosetta_scripts::ParsedProtocol peptide_parsed_protocol_mover;
	peptide_parsed_protocol_mover.parse_my_tag( mytag, data);                  // applying xml
	peptide_parsed_protocol_mover.apply(my_pose);                             // applying parsed protocol mover
	
}

// This function applies a user-defined fold tree if given
void apply_fold_tree_to_pose(std::string &fold_tree_file, std::string const &errmsg, core::pose::Pose &my_pose) {
	if (!utility::file::file_exists(fold_tree_file)) utility_exit_with_message(errmsg+"The file: "+fold_tree_file+" passed with the option -in:file:fold_tree does not exist in path.");
	utility::vector1< std::string > fold_tree_file_name_white_spacesplit_vector = utility::split_whitespace(fold_tree_file); // now check if this
	if (fold_tree_file_name_white_spacesplit_vector.size()==0) utility_exit_with_message(errmsg+"The application only found whitespaces given with the option -in:file:fold_tree so it will quit.");
	core::kinematics::FoldTreeOP usr_ftree(utility::pointer::make_shared< core::kinematics::FoldTree >());
	utility::io::izstream user_ftree_content_stream( fold_tree_file );
	user_ftree_content_stream >> *usr_ftree;
	user_ftree_content_stream.close();
	if (usr_ftree->check_fold_tree()){
		TR << "User Fold Tree is valid. Applying Fold Tree to Pose that was built from your XML tag" << std::endl;
		my_pose.fold_tree(*usr_ftree);
	} else utility_exit_with_message(errmsg+" User Fold Tree file is not valid. A good explanation of Fold Tree can be found in this link: https://www.rosettacommons.org/demos/latest/tutorials/fold_tree/fold_tree");
	
}

// This function reads the degrees of freedom file. The file will dictate how the built pose from the parsed protocol will be perturbed
void read_dof_file(std::string const &dofs_file_name, std::string const &errmsg, utility::vector1 < std::string >  &pertubation_all_array) {
	if (!utility::file::file_exists(dofs_file_name)) utility_exit_with_message(errmsg+"The file: "+dofs_file_name+" passed with the option -dofs_file does not exist in path.");
	utility::vector1< std::string > dofs_file_name_white_spacesplit_vector = utility::split_whitespace(dofs_file_name); // now check if this
	if (dofs_file_name_white_spacesplit_vector.size()==0) utility_exit_with_message(errmsg+"The application only found whitespaces given with the option -in:file:fold_tree so it will quit.");
	utility::io::izstream dofs_file_stream( dofs_file_name );
	utility::vector1 < std::string > dof_vector;
	std::string dof_line;
	while (getline(dofs_file_stream, dof_line)) pertubation_all_array.push_back(dof_line);
}

// This function reads the constraint file. The constraint file should have all the constraints that define the physical geometry of the system.
// Think about it as you are creating a custom ForceField for your pose.
ConstraintSetOP get_cst_set_from_cst_file(std::string const &cst_file, std::string const &errmsg, core::pose::Pose const &my_pose) {
	if (!utility::file::file_exists(cst_file)) utility_exit_with_message(errmsg+"The file: "+cst_file+" passed with the option -cst_file does not exist in path.");
	utility::vector1< std::string > cst_file_name_white_spacesplit_vector = utility::split_whitespace(cst_file); // now check if this
	if (cst_file_name_white_spacesplit_vector.size()==0) utility_exit_with_message(errmsg+"The application only found whitespaces given with the option -cst_file so it will quit.");
	std::string cst_file_content;
	utility::io::izstream cst_file_stream( cst_file);
	utility::slurp( cst_file_stream, cst_file_content );
	utility::vector1< std::string > cst_file_content_whitespace_split_vector = utility::split_whitespace(cst_file_content);
	if (cst_file_content_whitespace_split_vector.size()==0) utility_exit_with_message(errmsg+"Mesa your cst file is empty");
	return ConstraintIO::get_instance()->read_constraints( cst_file, utility::pointer::make_shared< ConstraintSet >(), my_pose ) ;
}

// This function will set up the qm score function. The score function that will be set up from this is going to be very dependent on the qm weights file that is being passed into.
// TODO: More descriptions on how to setup a qm score function needs to be added in Wiki
core::scoring::ScoreFunctionOP setup_qm_sfxn(std::string const & errmsg, std::string const & qm_weights_file) {
	if (!utility::file::file_exists(qm_weights_file)) utility_exit_with_message(errmsg+"The file: "+qm_weights_file+" passed with the option -score:weights does not exist in path.");
	utility::vector1< std::string > qm_weights_file_name_white_spacesplit_vector = utility::split_whitespace(qm_weights_file); // now check if this
	if (qm_weights_file_name_white_spacesplit_vector.size()==0) utility_exit_with_message(errmsg+"The application only found whitespaces given with the option -score:weights so it will quit.");
	return  core::scoring::ScoreFunctionFactory::create_score_function( qm_weights_file );
}

int
main( int argc, char * argv [] )
{
	try {
		register_options();
		devel::init( argc, argv );
		using namespace basic::thread_manager;
		// Start by instantiating a Pose Object. This an empty pose
		core::pose::Pose my_pose; 
		// 1. Build pose from the user-defined XML tag.
		std::string const xmltag_main_errmsg("Error in reading user defined XML tag: ");
		if (option[ pose_setup_file ].user()==false) utility_exit_with_message("A list of constraints for the crosslinker  is needed to build the pose, please try again by passing your constraints file with file option -cst_file");
		std::string xml_tag_file(option[ pose_setup_file ].value());
    build_initial_pose_from_xml_tag(xml_tag_file, xmltag_main_errmsg, my_pose);
		
		// 2. Apply user defined fold tree to pose built from the user-defined XML tag.
		std::string const ft_main_errmsg("Error in reading user defined Fold Tree: ");
    core::kinematics::FoldTreeOP pose_ftree(utility::pointer::make_shared< core::kinematics::FoldTree >());
		if (option[ in::file::fold_tree ].user()==false) TR.Info << "Since FoldTree is not passed by the user, the application will define the FoldTree from the pose." << std::endl;
        if (option[ in::file::fold_tree ].user()==true) {
            std::string fold_tree_file( option[ in::file::fold_tree ].value() );    
            apply_fold_tree_to_pose(fold_tree_file, ft_main_errmsg, my_pose);
        }
		
		// 3. Read in the dof file.
		std::string const dof_main_errmsg("Error in reading user defined Fold Tree: ");
		if (option[ dofs_file ].user()==false) utility_exit_with_message("User defined degrees of freedom (dof) is required for sampling the pose to be built, please try again by passing your dofs file with file option -dofs_file");
		std::string dofs_file_name( option[ dofs_file ].value() );
		utility::vector1 < std::string >   pertubation_all_array;
		read_dof_file(dofs_file_name, dof_main_errmsg, pertubation_all_array);
		
		// 4. Read in the cst file and return and setup contraint sfxn.
		std::string const cst_main_errmsg("Error in reading user defined constraint file: ");
		if (option[ constraints::cst_file ].user()==false) utility_exit_with_message(cst_main_errmsg + "User error. A list of constraints for the crosslinker  is needed to build the pose, please try again by passing your constraints file with file option -cst_file your_constraint_file");
		std::string cst_file( option[ constraints::cst_file ]()[1] ); // now cst_file_ is set to whatever the user has set it to.
		ConstraintSetOP cst_set ( get_cst_set_from_cst_file(cst_file, cst_main_errmsg, my_pose) );
		//core::scoring::ScoreFunctionOP constr_sfxn(setup_constr_sfxn());
		
		// 5. Read in the qm weights file and setup QM sfxn.
		std::string const qm_main_errmsg("Error in reading user defined QM weights file: ");
		if (option[ score::weights ].user()==false) utility_exit_with_message(qm_main_errmsg+"A score weights file with for QM weights is needed, please try again by passing your QM weights file with -score:weights");
		std::string qm_weights_file ( option[ score::weights ].value() );
		core::scoring::ScoreFunctionOP qm_sfxn(setup_qm_sfxn(qm_main_errmsg, qm_weights_file));
		
		// 6. Setup the sampling stage. 
		std::string const sampling_main_errmsg(" Error setting up the sampling stage: ");
		if (option[ pert_sample_number ].user()==false) utility_exit_with_message(sampling_main_errmsg+ "User must enter the total number of samples (integers) to be perturbed with the option -pert_sample_number");
		core::Size number_of_samples_for_perturbation( option[pert_sample_number].value() );
		if (number_of_samples_for_perturbation <= 0) utility_exit_with_message(sampling_main_errmsg+"Application exit in Sampling stage. User error. Sample size cannot be <= "+std::to_string(number_of_samples_for_perturbation));
		if (number_of_samples_for_perturbation == 1) TR.Warning << "Sampling size is small, however the protocol will continue. Consider setting it to at least over a 100";		
		PerturbationList pert_list_object(pertubation_all_array,my_pose);
	
   // 7. Setup fa_repulsion score function
   bool const user_has_provided_fa_rep_cutoff( option[ discard_sample_by_fa_repulsion_cutoff ].user() );
   core::Real const fa_repulsion_score_cuttoff( user_has_provided_fa_rep_cutoff ? option[ discard_sample_by_fa_repulsion_cutoff ].value() : 9999999999999999999999999999.9 );
   bool const user_has_provided_num_threads( option[ num_threads_qm ].user() );
   core::Size const number_of_threads_to_use( user_has_provided_num_threads ? option[ num_threads_qm ].value(): 0);
   //core::scoring::ScoreFunctionOP fa_rep_sfxn( setup_fa_rep_sfxn() );

    // 8. Do the sampling.
    bool do_dump_pdb = option[ dump ]();
    bool write_all_matrices = option [ debug_matrices ]();
		core::Real delta_qm_value = option [ delta_qm_energy_cutoff ]();
		PerturbationSampler pert_sampler_object(number_of_samples_for_perturbation,cst_set);
		pert_sampler_object.do_sampling(my_pose, pert_list_object, qm_sfxn, fa_repulsion_score_cuttoff, do_dump_pdb, write_all_matrices,delta_qm_value, number_of_threads_to_use);
		
    // 9. Do the linear algebra
    pert_sampler_object.lin_alg_least_square_fitting(write_all_matrices);

    // 10. Setup the user-defined output filename.
		std::string output_file_name(option[output_weights_file].value());

    // 11. Write report to user-defined filename and exit.
		pert_sampler_object.write_report_from_sampling(output_file_name);
 
    } catch ( utility::excn::Exception const & e ) {
    		e.display();
        return -1;
    	}
	return 0;
}
