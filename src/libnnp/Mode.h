// n2p2 - A neural network potential package
// Copyright (C) 2018 Andreas Singraber (University of Vienna)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef MODE_H
#define MODE_H

#include "CutoffFunction.h"
#include "Element.h"
#include "ElementMap.h"
#include "ErfcBuf.h"
#include "EwaldSetup.h"
#include "Log.h"
#include "ScreeningFunction.h"
#include "Settings.h"
#include "Structure.h"
#include "DFTD3.h"
#include "Repulsion.h"
#include "External.h"
#include "SymFnc.h"
#include <cstddef> // std::size_t
#include <map>     // std::map
#include <string>  // std::string
#include <vector>  // std::vector

namespace nnp
{

/** Base class for all NNP applications.
 *
 * This top-level class is the anchor point for existing and future
 * applications. It contains functions to set up an existing neural network
 * potential and calculate energies and forces for configurations given as
 * Structure. A minimal setup requires some consecutive functions calls as this
 * minimal example shows:
 *
 * ```
 * Mode mode;
 * mode.initialize();
 * mode.loadSettingsFile();
 * mode.setupElementMap();
 * mode.setupElements();
 * mode.setupCutoff();
 * mode.setupSymmetryFunctions();
 * mode.setupSymmetryFunctionGroups();
 * mode.setupSymmetryFunctionStatistics(false, false, true, false);
 * mode.setupNeuralNetwork();
 * ```
 * To load weights and scaling information from files add these lines:
 * ```
 * mode.setupSymmetryFunctionScaling();
 * mode.setupNeuralNetworkWeights();
 * ```
 * The NNP is now ready! If we load a structure from a data file:
 * ```
 * Structure structure;
 * ifstream file;
 * file.open("input.data");
 * structure.setupElementMap(mode.elementMap);
 * structure.readFromFile(file);
 * file.close();
 * ```
 * we can finally predict the energy and forces from the neural network
 * potential:
 * ```
 * structure.calculateNeighborList(mode.getMaxCutoffRadius());
 * mode.calculateSymmetryFunctionGroups(structure, true);
 * mode.calculateAtomicNeuralNetworks(structure, true);
 * mode.calculateEnergy(structure);
 * mode.calculateForces(structure);
 * cout << structure.energy << '\n';
 * ```
 * The resulting potential energy is stored in Structure::energy, the forces
 * on individual atoms are located within the Structure::atoms vector in
 * Atom::f.
 */
class Mode
{
public:
    enum class NNPType
    {
        /// Short range NNP (2G-HDNNP).
        HDNNP_2G = 2,
        /// NNP with electrostatics and non-local charge transfer (4G-HDNNP).
        HDNNP_4G = 4,
        /** Short range NNP with charge NN, no electrostatics/Qeq (M. Bircher).
         *
         * This is a simplified version of a 4G-HDNNP. Two neural networks are
         * used: the first one predicts atomic charges, which will be used for
         * the second NN as additional input neuron. There is no electrostatic
         * energy and no global charge equilibration as in 4G-HDNNP.
         */
        HDNNP_Q = 10
    };

    Mode();
    /** Write welcome message with version information.
     */
    void                     initialize();
    /** Open settings file and load all keywords into memory.
     *
     * @param[in] fileName Settings file name.
     */
    void                     loadSettingsFile(std::string const& fileName
                                                                 = "input.nn");
    /** Combine multiple setup routines and provide a basic NNP setup.
     *
     * @param[in] nnpDir Optional directory where NNP files reside.
     * @param[in] initialHardness Signalizes to use initial hardness in NN
     *                  settings.
     *
     * Sets up elements, symmetry functions, symmetry function groups, neural
     * networks. No symmetry function scaling data is read, no weights are set.
     */
    void                     setupGeneric(std::string const& nnpDir = "",
                                          bool      initialHardness = false);
    /** Set up global.
     *
     * If the keywords `verbose` is present, we print more informations
     */
    void                     setupGlobal();

    /** Set up external.
     *
     * If the keywords `external_progname` is present, we add external energy (&forces)
     * User can set type prefixname of input/ouput files of external_progname) using external_prefixfilename
     */
    void                     setupExternal();

    /** Set up Repulsion.
     *
     * If the keywords `repulsion_type` is present, we add repulsion energy (&forces)
     * User can set type o potenrion using repulsion tyepRepulsion(ZBL, Moliere, KrC)
     * User can set rIn using repulsion_rin keyword (default rIn=0)
     * User can set rOut using repulsion_rout keyword(default rOut=rcut)
     */
    void                     setupRepulsion();
    /** Set up DFT-D3.
     *
     * If the keywords `dftd3_func` is present, we add DFTD3 energy (&forces)
     * User can set Ithe version using dftd3_version keyword
     * User can set Ithe version using dftd3_verbose keyword
     * We ASSUME THAT ALL DATA are in atomic unit
     */
    void                     setupDFTD3();
    /** Set up normalization.
     *
     * If the keywords `mean_energy`, `conv_length` and
     * `conv_length` are present, the provided conversion factors are used to
     * internally use a different unit system.
     */
    void                     setupNormalization();
    /** Set up the element map.
     *
     * Uses keyword `elements`. This function should follow immediately after
     * settings are loaded via loadSettingsFile().
     */
    virtual void             setupElementMap();
    /** Set up all Element instances.
     *
     * Uses keywords `number_of_elements` and `atom_energy`. This function
     * should follow immediately after setupElementMap().
     */
    virtual void             setupElements();
    /** Set up cutoff function for all symmetry functions.
     *
     * Uses keyword `cutoff_type`. Cutoff parameters are read from settings
     * keywords and stored internally. As soon as setupSymmetryFunctions() is
     * called the settings are restored and used for all symmetry functions.
     * Thus, this function must be called before setupSymmetryFunctions().
     */
    void                     setupCutoff();
    /** Set up all symmetry functions.
     *
     * Uses keyword `symfunction_short`. Reads all symmetry functions from
     * settings and automatically assigns them to the correct element.
     */
    virtual void             setupSymmetryFunctions();
    /** Set up "empy" symmetry function scaling.
     *
     * Does not use any keywords. Sets no scaling for all symmetry functions.
     * Call after setupSymmetryFunctions(). Alternatively set scaling via
     * setupSymmetryFunctionScaling().
     */
    void                     setupSymmetryFunctionScalingNone();
    /** Set up symmetry function scaling from file.
     *
     * @param[in] fileName Scaling file name.
     *
     * Uses keywords `scale_symmetry_functions`, `center_symmetry_functions`,
     * `scale_symmetry_functions_sigma`, `scale_min_short` and
     * `scale_max_short`. Reads in scaling information and sets correct scaling
     * behavior for all symmetry functions. Call after
     * setupSymmetryFunctions().
     */
    virtual void             setupSymmetryFunctionScaling(
                                 std::string const& fileName = "scaling.data");
    /** Set up symmetry function groups.
     *
     * Does not use any keywords. Call after setupSymmetryFunctions() and
     * ensure that correct scaling behavior has already been set.
     */
    virtual void             setupSymmetryFunctionGroups();
#ifndef NNP_NO_SF_CACHE
    /** Set up symmetry function cache.
     *
     * @param[in] verbose If true, print more cache information.
     *
     * Searches symmetry functions for identical cutoff functions or compact
     * function (i.e. all cachable stuff) and sets up a caching index.
     */
    virtual void             setupSymmetryFunctionCache(bool verbose = false);
#endif
    /** Extract required memory dimensions for symmetry function derivatives.
     *
     * @param[in] verbose If true, print all symmetry function lines.
     *
     * Call after symmetry functions have been set up and sorted.
     */
    void                     setupSymmetryFunctionMemory(bool verbose = false);
    /** Set up symmetry function statistics collection.
     *
     * @param[in] collectStatistics Whether statistics (min, max, mean, sigma)
     *                              is collected.
     * @param[in] collectExtrapolationWarnings Whether extrapolation warnings
     *                                         are logged.
     * @param[in] writeExtrapolationWarnings Write extrapolation warnings
     *                                       immediately when they occur.
     * @param[in] stopOnExtrapolationWarnings Throw error immediately when
     *                                        an extrapolation warning occurs.
     *
     * Does not use any keywords. Calling this setup function is not required,
     * by default no statistics collection is enabled (all arguments `false`).
     * Call after setupElements().
     */
    void                     setupSymmetryFunctionStatistics(
                                             bool collectStatistics,
                                             bool collectExtrapolationWarnings,
                                             bool writeExtrapolationWarnings,
                                             bool stopOnExtrapolationWarnings);
    /** Setup matrix storing all symmetry function cut-offs for each element.
     */
    void                     setupCutoffMatrix();
    /** Set up neural networks for all elements.
     *
     * Uses keywords `global_hidden_layers_short`, `global_nodes_short`,
     * `global_activation_short`, `normalize_nodes`. Call after
     * setupSymmetryFunctions(), only then the number of input layer neurons is
     * known.
     */
    virtual void             setupNeuralNetwork();
    /** Set up neural network weights from files with given name format.
     *
     * @param[in] fileNameFormats Map of NN ids to format for weight file
     *                            names. Must contain a placeholder "%zu" for
     *                            the element atomic number. Map keys are
     *                            "short", "elec". Map argument may be
     *                            ommitted, then default name formats are used.
     *
     * Does not use any keywords. The weight files should contain one weight
     * per line, see NeuralNetwork::setConnections() for the correct order.
     */
    virtual void             setupNeuralNetworkWeights(
                                 std::map<std::string,
                                          std::string> fileNameFormats =
                                 std::map<std::string,
                                          std::string>());
    /** Set up neural network weights from files with given name format.
     *
     * @param[in] directoryPrefix Directory prefix which is applied to all
     *                            fileNameFormats.
     * @param[in] fileNameFormats Map of NN ids to format for weight file
     *                            names. Must contain a placeholder "%zu" for
     *                            the element atomic number. Map keys are
     *                            "short", "elec". Map argument may be
     *                            ommitted, then default name formats are used.
     *
     * Does not use any keywords. The weight files should contain one weight
     * per line, see NeuralNetwork::setConnections() for the correct order.
     */
    virtual void             setupNeuralNetworkWeights(
                                 std::string           directoryPrefix,
                                 std::map<std::string,
                                          std::string> fileNameFormats =
                                 std::map<std::string,
                                          std::string>());
    /** Set up electrostatics related stuff (hardness, screening, ...).
     *
     * @param[in] initialHardness Use initial hardness from keyword in settings
     *                            file (useful for training).
     * @param[in] directoryPrefix Directory prefix which is applied to
     *                            fileNameFormat.
     * @param[in] fileNameFormat Name format of file containing atomic
     *                           hardness data.
     */
    virtual void             setupElectrostatics(bool initialHardness =
                                                 false,
                                                 std::string directoryPrefix =
                                                 "",
                                                 std::string fileNameFormat =
                                                 "hardness.%03zu.data");
    /** Calculate all symmetry functions for all atoms in given structure.
     *
     * @param[in] structure Input structure.
     * @param[in] derivatives If `true` calculate also derivatives of symmetry
     *                        functions.
     *
     * This function should be replaced by calculateSymmetryFunctionGroups()
     * whenever possible. Results are stored in Atom::G. If derivatives are
     * calculated, additional results are stored in Atom::dGdr and
     * Atom::Neighbor::dGdr.
     */
    void                     calculateSymmetryFunctions(
                                                       Structure& structure,
                                                       bool const derivatives);
    /** Calculate all symmetry function groups for all atoms in given
     * structure.
     *
     * @param[in] structure Input structure.
     * @param[in] derivatives If `true` calculate also derivatives of symmetry
     *                        functions.
     *
     * This function replaces calculateSymmetryFunctions() when symmetry
     * function groups are enabled (faster, default behavior). Results are
     * stored in Atom::G. If derivatives are calculated, additional results are
     * stored in Atom::dGdr and Atom::Neighbor::dGdr.
     */
    void                     calculateSymmetryFunctionGroups(
                                                       Structure& structure,
                                                       bool const derivatives);
    // /** Calculate a single atomic neural network for a given atom and nn type.
    // *
    // * @param[in] nnId Neural network identifier, e.g. "short", "charge".
    // * @param[in] atom Input atom.
    // * @param[in] derivatives If `true` calculate also derivatives of neural
    // *                        networks with respect to input layer neurons
    // *                        (required for force calculation).
    // *
    // * The atomic energy and charge is stored in Atom::energy and Atom::charge,
    // * respectively. If derivatives are calculated the results are stored in
    // * Atom::dEdG or Atom::dQdG.
    // */
    //void                     calculateAtomicNeuralNetwork(
    //                                           std::string const& nnId,
    //                                           Atom&              atom,
    //                                           bool const         derivatives);
    /** Calculate atomic neural networks for all atoms in given structure.
     *
     * @param[in] structure Input structure.
     * @param[in] derivatives If `true` calculate also derivatives of neural
     *                        networks with respect to input layer neurons
     *                        (required for force calculation).
     * @param[in] id Neural network ID to use. If empty, the first entry
     *                        nnk.front() is used.
     * @param[in] suppressOutput Turn on/off output generation (for example
     *                        during training.
     */
    void                     calculateAtomicNeuralNetworks(
                                           Structure&  structure,
                                           bool const  derivatives,
                                           std::string id = "",
                                           bool const  suppressOutput = false);
    /** Perform global charge equilibration method.
     *
     * @param[in] structure Input structure.
     * @param[in] derivativesElec Turn on/off calculation of dElecdQ and
     *                          pElecpr (Typically needed for elecstrosttic
     *                          forces).
     * @param[in] suppressOutput Turn on/off output generation (for example
     *                          during training.
     */
    void                     chargeEquilibration(
                                                Structure& structure,
                                                bool const derivativesElec,
                                                bool const suppressOutput = false);
    /** Calculate potential energy for a given structure.
     *
     * @param[in] structure Input structure.
     *
     * Sum up potential energy from atomic energy contributions. Result is
     * stored in Structure::energy.
     */
    void                     calculateEnergy(Structure& structure) const;
    /** Calculate total charge for a given structure.
     *
     * @param[in] structure Input structure.
     *
     * Sum up charge from atomic charge contributions. Result is
     * stored in Structure::charge.
     */
    void                     calculateCharge(Structure& structure) const;
    /** Calculate forces for all atoms in given structure.
     *
     * @param[in] structure Input structure.
     *
     * Combine intermediate results from symmetry function and neural network
     * computation to atomic forces. Results are stored in Atom::f.
     */
    void                     calculateForces(Structure& structure) const;
    /** Evaluate neural network potential (includes total energy, optionally
     *  forces and in some cases charges.
     *  @param[in] structure Input structure.
     *  @param[in] useForces If true, calculate forces too.
     *  @param[in] useDEdG If true, calculate dE/dG too.
     */
    void                     evaluateNNP(Structure& structure,
                                         bool useForces = true,
                                         bool useDEdG = true);
    /** Add atomic energy offsets to reference energy.
     *
     * @param[in] structure Input structure.
     * @param[in] ref If true, use reference energy, otherwise use NN energy.
     */
    void                     addEnergyOffset(Structure& structure,
                                             bool       ref = true);
    /** Remove atomic energy offsets from reference energy.
     *
     * @param[in] structure Input structure.
     * @param[in] ref If true, use reference energy, otherwise use NN energy.
     *
     * This function should be called immediately after structures are read in.
     */
    void                     removeEnergyOffset(Structure& structure,
                                                bool       ref = true);
    /** Get atomic energy offset for given structure.
     *
     * @param[in] structure Input structure.
     *
     * @return Summed atomic energy offsets for structure.
     */
    double                   getEnergyOffset(Structure const& structure) const;
    /** Add atomic energy offsets and return energy.
     *
     * @param[in] structure Input structure.
     * @param[in] ref If true, use reference energy, otherwise use NN energy.
     *
     * @return Reference or NNP energy with energy offsets added.
     *
     * @note If normalization is used, ensure that structure energy is already
     * in physical units.
     */
    double                   getEnergyWithOffset(
                                            Structure const& structure,
                                            bool             ref = true) const;
    /** Apply normalization to given property.
     *
     * @param[in] property One of "energy", "force", "charge", "hardness"
     *                     "negativity".
     * @param[in] value Input property value in physical units.
     *
     * @return Property in normalized units.
     */
    double                   normalized(std::string const& property,
                                        double             value) const;
    /** Apply normalization to given energy of structure.
     *
     * @param[in] structure Input structure with energy in physical units.
     * @param[in] ref If true, use reference energy, otherwise use NN energy.
     *
     * @return Energy in normalized units.
     */
    double                   normalizedEnergy(
                                            Structure const& structure,
                                            bool             ref = true) const;
    /** Undo normalization for a given property.
     *
     * @param[in] property One of "energy", "force", "charge", "hardness",
     *                     "negativity".
     * @param[in] value Input property value in normalized units.
     *
     * @return Property in physical units.
     */
    double                   physical(std::string const& property,
                                      double             value) const;
    /** Undo normalization for a given energy of structure.
     *
     * @param[in] structure Input structure with energy in normalized units.
     * @param[in] ref If true, use reference energy, otherwise use NN energy.
     *
     * @return Energy in physical units.
     */
    double                   physicalEnergy(Structure const& structure,
                                            bool             ref = true) const;
    /** Convert one structure to normalized units.
     *
     * @param[in,out] structure Input structure.
     */
    void                     convertToNormalizedUnits(
                                                   Structure& structure) const;
    /** Convert one structure to physical units.
     *
     * @param[in,out] structure Input structure.
     */
    void                     convertToPhysicalUnits(
                                                   Structure& structure) const;
    /** Logs Ewald params whenever they change.
     *
     */
    void                     logEwaldCutoffs();
    /** Count total number of extrapolation warnings encountered for all
     * elements and symmetry functions.
     *
     * @return Number of extrapolation warnings.
     */
    std::size_t              getNumExtrapolationWarnings() const;
    /** Erase all extrapolation warnings and reset counters.
     */
    void                     resetExtrapolationWarnings();
    /** Getter for Mode::nnpType.
     *
     * @return HDNNP type (2G, 4G,..)
     */
    NNPType                  getNnpType() const;
    /** Getter for Mode::meanEnergy.
     *
     * @return Mean energy per atom.
     */
    double                   getMeanEnergy() const;
    /** Getter for Mode::convEnergy.
     *
     * @return Energy unit conversion factor.
     */
    double                   getConvEnergy() const;
    /** Getter for Mode::convLength.
     *
     * @return Length unit conversion factor.
     */
    double                   getConvLength() const;
    /** Getter for Mode::convCharge.
     *
     * @return Charge unit conversion factor.
     */
    double                   getConvCharge() const;
    /** Getter for Mode::convEnergyToHartree.
     *
     * @return Energy unit conversion factor to Hartree (needed for D3, Rep and External).
     */
    double                   getConvToHartree() const;
    /** Getter for Mode::convLengthToBohr.
     *
     * @return Length unit conversion factor to Bohr (needed for D3, Rep and External).
     */
    double                   getConvToBohr() const;
    /** Getter for Mode::maxCutoffRadius.
     *
     * @return Maximum cutoff radius of all symmetry functions.
     *
     * The maximum cutoff radius is determined by setupSymmetryFunctions().
     */
    double                   getMaxCutoffRadius() const;
    /** Getter for Mode::numElements.
     *
     * @return Number of elements defined.
     *
     * The number of elements is determined by setupElements().
     */
    std::size_t              getNumElements() const;
    /** Get number of symmetry functions per element.
     *
     * @return Vector with number of symmetry functions for each element.
     */
    std::vector<std::size_t> getNumSymmetryFunctions() const;
    /** Check if normalization is enabled.
     *
     * @return Value of #normalize.
     */
    bool                     useNormalization() const;
    /** Check if keyword was found in settings file.
     *
     * @param[in] keyword Keyword for which value is requested.
     *
     * @return `true` if keyword exists, `false` otherwise.
     */
    bool                     settingsKeywordExists(
                                             std::string const& keyword) const;
    /** Get value for given keyword in Settings instance.
     *
     * @param[in] keyword Keyword for which value is requested.
     *
     * @return Value string corresponding to keyword.
     */
    std::string              settingsGetValue(
                                             std::string const& keyword) const;
    /** Prune symmetry functions according to their range and write settings
     * file.
     *
     * @param[in] threshold Symmetry functions with range (max - min) smaller
     *                      than this threshold will be pruned.
     *
     * @return List of line numbers with symmetry function to be removed.
     */
    std::vector<std::size_t> pruneSymmetryFunctionsRange(double threshold);
    /** Prune symmetry functions with sensitivity analysis data.
     *
     * @param[in] threshold Symmetry functions with sensitivity lower than this
     *                      threshold will be pruned.
     * @param[in] sensitivity Sensitivity data for each element and symmetry
     *                        function.
     *
     * @return List of line numbers with symmetry function to be removed.
     */
    std::vector<std::size_t> pruneSymmetryFunctionsSensitivity(
                                            double threshold,
                                            std::vector<
                                            std::vector<double> > sensitivity);
    /** Copy settings file but comment out lines provided.
     *
     * @param[in] prune List of line numbers to comment out.
     * @param[in] fileName Output file name.
     */
    void                     writePrunedSettingsFile(
                                              std::vector<std::size_t> prune,
                                              std::string              fileName
                                                          = "output.nn") const;
    /** Write complete settings file.
     *
     * @param[in,out] file Settings file.
     */
    void                     writeSettingsFile(
                                             std::ofstream* const& file) const;

    /// Global element map, populated by setupElementMap().
    ElementMap elementMap;
    /// Global log file.
    Log        log;

protected:
    /// Setup data for one neural network.
    struct NNSetup
    {
        struct Topology
        {
            /// Number of NN layers (including input and output layer).
            int                                numLayers;
            /// Number of neurons per layer.
            std::vector<int>                   numNeuronsPerLayer;
            /// Activation function type per layer.
            std::vector<
            NeuralNetwork::ActivationFunction> activationFunctionsPerLayer;

            /// Constructor.
            Topology() : numLayers(0) {};
        };

        /// NN identifier, e.g. "short", "charge",...
        std::string           id;
        /// Description string for log output, e.g. "electronegativity".
        std::string           name;
        /// Format for weight files.
        std::string           weightFileFormat;
        /// Suffix for keywords (NN topology related).
        std::string           keywordSuffix;
        /// Suffix for some other keywords (weight file loading related).
        std::string           keywordSuffix2;
        /// Per-element NN topology.
        std::vector<Topology> topology;
    };


    NNPType                    nnpType;
    bool                       normalize;
    bool                       checkExtrapolationWarnings;
    std::size_t                numElements;
    std::vector<std::size_t>   minNeighbors;
    std::vector<double>        minCutoffRadius;
    double                     maxCutoffRadius;
    double                     cutoffAlpha;
    double                     meanEnergy;
    double                     convEnergy;
    double                     convLength;
    double                     convCharge;
    double                     convToHartree;
    double                     convToBohr;
    int                        verbose;
    EwaldSetup                 ewaldSetup;
    settings::Settings         settings;
    SymFnc::ScalingType        scalingType;
    CutoffFunction::CutoffType cutoffType;
    ScreeningFunction          screeningFunction;
    std::vector<Element>       elements;
    std::vector<std::string>   nnk;
    std::map<
    std::string, NNSetup>      nns;
    /// Matrix storing all symmetry function cut-offs for all elements.
    std::vector<
    std::vector<double>>       cutoffs;
    ErfcBuf                    erfcBuf;
    DFTD3                    dftD3;
    Repulsion                rep;
    External                 external;

    /** Read in weights for a specific type of neural network.
     *
     * @param[in] id Actual network type to initialize ("short" or "elec").
     * @param[in] fileNameFormat Weights file name format.
     */
    void readNeuralNetworkWeights(std::string const& id,
                                  std::string const& fileName);
};

//////////////////////////////////
// Inlined function definitions //
//////////////////////////////////

inline Mode::NNPType Mode::getNnpType() const
{
    return nnpType;
}

inline double Mode::getMeanEnergy() const
{
    return meanEnergy;
}

inline double Mode::getConvEnergy() const
{
    return convEnergy;
}

inline double Mode::getConvLength() const
{
    return convLength;
}

inline double Mode::getConvCharge() const
{
    return convCharge;
}
inline double Mode::getConvToHartree() const
{
    return convToHartree;
}

inline double Mode::getConvToBohr() const
{
    return convToBohr;
}

inline double Mode::getMaxCutoffRadius() const
{
    return maxCutoffRadius;
}

inline std::size_t Mode::getNumElements() const
{
    return numElements;
}

inline bool Mode::useNormalization() const
{
    return normalize;
}

}

#endif
