// This file is part of the ExaHyPE2 project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#pragma once


#include "exahype2/CellData.h"
#include "exahype2/enumerator/AoSLexicographicEnumerator.h"


namespace exagrype {
  namespace ccz4 {
      /**
       * Recompute auxiliary variables for FD4 scheme with a 4th order scheme
       *
       * We are given a patch as well as the number of mesh cells per axis
       * within this patch. I also parameterise over the unknowns and the
       * auxiliaryVariables. By default, we'd expect 58 and 0 here, but we
       * want to routine to work for augmented systems, too.
       *
       * The reconstruction uses a 4th order Finite Difference scheme. The
       * implementation in exagrype::ccz4::internal::recomputeAuxiliaryVariablesFD4_4thOrder_LoopBody()
       * realises the normal stencil
       *
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       * [1 -8 0 8 -1]
       * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       *
       * scaled with @f$ \frac{1}{12h} @f$. The stencil reaches two layers
       * around every grid cell, so the routine computes the derivatives of the
       * patch's own grid cells only and needs a halo of two.
       *
       * ## Why not one layer outside of the patch
       *
       * @image html recomputeAuxiliaryVariablesFD4_4thOrder.png
       *
       * An earlier version also evaluated the derivatives one layer outside of
       * the patch (brownish layer to the right), with the argument that the
       * fluxes and non-conservative products on the patch's faces need the
       * derivatives on both sides. This required a halo of three. The Finite
       * Differences of exahype2.api.rkfd4 evaluate no terms on faces, though:
       * The derivatives of the halo are computed by the neighbouring patch,
       * which owns these grid cells, and are exchanged via the faces
       * (exagrype.api.FOSOReconstructionSolver). Furthermore, Peano holds no
       * diagonal neighbours: What we really get is not the illustration above
       * to the left but the one in the centre. The diagonal blocks hold
       * garbage, so in the outer layer only the derivatives along its normal
       * would be correct.
       *
       *
       * @param patchData          Host the actual patch data. For this
       *   function, we only have to ensure that the QIn (data plus halo)
       *   are properly set, and that the patch size in patchData is
       *   correct. All the other (meta data) properties have no influence
       *   and can hold any data.
       * @param numberOfGridCellsPerPatchPerAxis Has to be bigger or equal to
       *   haloSize.
       * @param haloSize           Has to be at least 2.
       * @param unknowns           Typically 58, but can be bigger if you hold
       *   additional quantities within the PDE.
       * @param auxiliaryVariables Typically 0.
       *
       * @see exagrype::ccz4::internal::recomputeAuxiliaryVariablesFD4_4thOrder_LoopBody()
       */
      void recomputeAuxiliaryVariablesFD4_4thOrder(
        ::exahype2::CellData<double, double>&   patchData,
        int                                     numberOfGridCellsPerPatchPerAxis,
        int                                     haloSize,
        int                                     unknowns,
        int                                     auxiliaryVariables
      );

      void recomputeAuxiliaryVariablesFD4_centralDifferences(
        ::exahype2::CellData<double, double>&   patchData,
        int                                     numberOfGridCellsPerPatchPerAxis,
        int                                     haloSize,
        int                                     unknowns,
        int                                     auxiliaryVariables
      );

      void recomputeAuxiliaryVariablesFD4_leftDifferences(
        ::exahype2::CellData<double, double>&   patchData,
        int                                     numberOfGridCellsPerPatchPerAxis,
        int                                     haloSize,
        int                                     unknowns,
        int                                     auxiliaryVariables
      );

      void recomputeAuxiliaryVariablesFD4_rightDifferences(
        ::exahype2::CellData<double, double>&   patchData,
        int                                     numberOfGridCellsPerPatchPerAxis,
        int                                     haloSize,
        int                                     unknowns,
        int                                     auxiliaryVariables
      );

      /**
       * Reconstruct the 33 first derivatives of the 25 evolved CCZ4 unknowns
       *
       * This is the entry point of exagrype.api.FOSOReconstructionSolver. The
       * second-order solver evolves the first 25 of the 58 CCZ4 entries
       * (conformal metric, extrinsic curvature, ..., phi); the remaining 33
       * entries are the first derivatives of the lapse, shift, metric and
       * phi. This routine computes these derivatives (entries 25 to 57) of
       * one cell from the 25 evolved unknowns with their halo.
       *
       * The routine is a thin wrapper around
       * recomputeAuxiliaryVariablesFD4_4thOrder(), i.e. the fourth-order
       * stencil @f$ [1\ -8\ 0\ 8\ -1]/(12h) @f$. That routine expects the
       * 58 entries per grid cell interleaved, and a halo of at least two.
       * The wrapper assembles such a buffer, calls the routine and copies
       * the derivatives of the cell's own grid cells back.
       *
       * The routine computes the derivatives of the cell's own grid cells
       * only. The CCZ4 solver obtains the derivatives of its halo from the
       * neighbouring cells via the faces, as for the Finite Volumes
       * (reconstructFirstDerivativesCentralDifferences()). A halo of two, the
       * one of the plain Finite Differences, hence suffices.
       *
       * @param primaryVariablesWithHalo The 25 evolved unknowns per grid
       *   cell of the cell with its halo, i.e.
       *   (numberOfGridCellsPerPatchPerAxis+2*haloSize)^3 grid cells in the
       *   layout of a reconstructed patch (peano4::utils::dLinearised()).
       * @param patchSize Extent of the cell (without halo).
       * @param numberOfGridCellsPerPatchPerAxis Grid cells per axis within
       *   the cell.
       * @param haloSize Has to be at least two.
       * @param firstDerivatives Output, 33 entries per grid cell of the cell
       *   (without halo) in the layout of a patch. Entry k corresponds to
       *   entry 25+k of the 58 CCZ4 unknowns.
       */
      void reconstructFirstDerivativesFD4(
        const double* __restrict__                  primaryVariablesWithHalo,
        const tarch::la::Vector<Dimensions,double>& patchSize,
        int                                         numberOfGridCellsPerPatchPerAxis,
        int                                         haloSize,
        double* __restrict__                        firstDerivatives
      );

      /**
       * Reconstruct the 33 first derivatives with second-order central differences
       *
       * Same contract as reconstructFirstDerivativesFD4(), but the stencil is
       * @f$ [-1\ 0\ 1]/(2h) @f$ (recomputeAuxiliaryVariablesFD4_centralDifferences()).
       * It reaches one layer out only, so a halo of one grid cell suffices.
       * This is what the Finite Volumes (overlap one) of
       * exagrype.api.FOSOCCZ4Solver use: the derivatives of the halo layer are
       * not computed here but exchanged via the faces, as they are computed by the
       * neighbouring cell.
       *
       * @param primaryVariablesWithHalo See reconstructFirstDerivativesFD4().
       * @param patchSize See reconstructFirstDerivativesFD4().
       * @param numberOfGridCellsPerPatchPerAxis See reconstructFirstDerivativesFD4().
       * @param haloSize Has to be at least one.
       * @param firstDerivatives See reconstructFirstDerivativesFD4().
       */
      void reconstructFirstDerivativesCentralDifferences(
        const double* __restrict__                  primaryVariablesWithHalo,
        const tarch::la::Vector<Dimensions,double>& patchSize,
        int                                         numberOfGridCellsPerPatchPerAxis,
        int                                         haloSize,
        double* __restrict__                        firstDerivatives
      );

      namespace internal {
        /**
         * Recompute auxiliary variables
         *
         * This function calculates the auxiliary variables as gradietns of the primary
         * variables. The graident is calculated using a 5-point stencil in 1D which has
         * 2 neighbouring cells on each side of the central one. The computation requires
         * us to have access to 2 layers around each grid cell. The caller
         * (recomputeAuxiliaryVariablesFD4_4thOrder()) evaluates it on the
         * patch's own grid cells, i.e. it needs a halo of two.
         *
         */
        void recomputeAuxiliaryVariablesFD4_4thOrder_LoopBody(
          double* __restrict__                                       QIn,
          const ::exahype2::enumerator::AoSLexicographicEnumerator& QInEnumerator,
          const tarch::la::Vector<Dimensions,double>&  patchCentre,
          const tarch::la::Vector<Dimensions,double>&  patchSize,
          int                                          patchIndex,
          const tarch::la::Vector<Dimensions,int>&     volumeIndex,
          int                                          normal
        );


        /**
         *
         */
        void recomputeAuxiliaryVariablesFD4_centralDifferences_LoopBody(
          double* __restrict__                                       QIn,
          const ::exahype2::enumerator::AoSLexicographicEnumerator& QInEnumerator,
          const tarch::la::Vector<Dimensions,double>&  patchCentre,
          const tarch::la::Vector<Dimensions,double>&  patchSize,
          int                                          patchIndex,
          const tarch::la::Vector<Dimensions,int>&     volumeIndex,
          int                                          normal
        );


        /**
         *
         */
        void recomputeAuxiliaryVariablesFD4_leftDifferences_LoopBody(
          double* __restrict__                                       QIn,
          const ::exahype2::enumerator::AoSLexicographicEnumerator& QInEnumerator,
          const tarch::la::Vector<Dimensions,double>&  patchCentre,
          const tarch::la::Vector<Dimensions,double>&  patchSize,
          int                                          patchIndex,
          const tarch::la::Vector<Dimensions,int>&     volumeIndex,
          int                                          normal
        );


        /**
         *
         */
        void recomputeAuxiliaryVariablesFD4_rightDifferences_LoopBody(
          double* __restrict__                                       QIn,
          const ::exahype2::enumerator::AoSLexicographicEnumerator& QInEnumerator,
          const tarch::la::Vector<Dimensions,double>&  patchCentre,
          const tarch::la::Vector<Dimensions,double>&  patchSize,
          int                                          patchIndex,
          const tarch::la::Vector<Dimensions,int>&     volumeIndex,
          int                                          normal
        );
      }
  }
}
