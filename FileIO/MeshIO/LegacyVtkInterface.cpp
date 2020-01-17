/**
 * \file LegacyVtkInterface.cpp
 * 05/04/2011 LB Initial implementation
 *
 * Implementation of LegacyVtkInterface class
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ** INCLUDES **
#include "LegacyVtkInterface.h"

#include "FEMEnums.h"
#include "ProcessInfo.h"
#include "fem_ele_std.h"
#include "matrix_class.h"
#include "msh_lib.h"
#include "msh_mesh.h"
#include "rf_mmp_new.h"  // this is for class CMediumProperties, what else???
#include "rf_pcs.h"
#include "rf_pcs.h"

#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif  // GEM_REACT

#include <iomanip>
#include <string>

#include "FileTools.h"
#include "Output.h"

using namespace std;

LegacyVtkInterface::LegacyVtkInterface(MeshLib::CFEMesh* mesh,
                                       std::vector<std::string>
                                           pointArrayNames,
                                       std::vector<std::string>
                                           cellArrayNames,
                                       std::vector<std::string>
                                           materialPropertyArrayNames,
                                       std::string meshTypeName,
                                       ProcessInfo* processInfo)
    : _mesh(mesh),
      _pointArrayNames(pointArrayNames),
      _cellArrayNames(cellArrayNames),
      _materialPropertyArrayNames(materialPropertyArrayNames),
      _meshTypeName(meshTypeName),
      _processInfo(processInfo)
{
    _processType = convertProcessTypeToString(processInfo->getProcessType());
    _mesh = FEMGet(_processType);
}

LegacyVtkInterface::~LegacyVtkInterface() {}

#if defined(USE_PETSC)
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2006 KG44 Implementation
   09/2006 KG44 Output for MPI - correct OUTPUT not yet implemented
   12/2008 NW Remove ios::app, Add PCS name to VTK file name
**************************************************************************/
void LegacyVtkInterface::WriteDataVTKPETSC(int number,
                                           double simulation_time,
                                           std::string baseFilename) const
{
    if (!_mesh)
    {
        cout << "Warning in LegacyVtkInterface::WriteVTKNodes - no MSH data"
             << "\n";
        return;
    }

    if (_processType.compare("INVALID_PROCESS") != 0)
        baseFilename += "_" + _processType;

    stringstream ss;
    // setw(4) sets the number of digits to be used
    // and creates leading zeros if necessary
    ss << setw(4) << setfill('0') << number;
    baseFilename += ss.str();
    baseFilename += ".vtk";

    // convert filename string to *char
    const char* filename = baseFilename.c_str();

    // Open the viewer
    PetscViewer viewer;
    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII);
    PetscViewerFileSetName(viewer, filename);
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
    // write header
    PetscViewerASCIIPrintf(viewer, "# vtk DataFile Version 3.0\n");
    PetscViewerASCIIPrintf(viewer, "Unstructured Grid from OpenGeoSys-GEM\n");
    PetscViewerASCIIPrintf(viewer, "ASCII\n");
    PetscViewerASCIIPrintf(viewer, "DATASET UNSTRUCTURED_GRID\n");
    PetscViewerASCIIPrintf(viewer, "FIELD TimesAndCycles 2\n");
    PetscViewerASCIIPrintf(viewer, "TIME 1 1 double\n");
    PetscViewerASCIIPrintf(viewer, "%lg\n", simulation_time);
    PetscViewerASCIIPrintf(viewer, "CYCLE 1 1 long\n");
    PetscViewerASCIIPrintf(viewer, "%d\n", number);

    this->WriteVTKPointDataPETSC(viewer);
    this->WriteVTKCellDataPETSC(viewer);
    this->WriteVTKDataArraysPETSC(viewer);

    // close viewer
    PetscViewerDestroy(&viewer);

    return;
}

void LegacyVtkInterface::WriteVTKPointDataPETSC(PetscViewer viewer) const
{
    PetscScalar *xp, *yp, *zp;  // used for pointer
    PetscInt low, high, nn;
    PetscInt count;
    int i;
    VecScatter ctx;

    MeshLib::CFEMesh* mesh = fem_msh_vector[0];
    //    const int nn = mesh->getNumNodesGlobal(); //global number of nodes
    //    without shadow nodes
    const size_t n_linear_pnts(mesh->GetNodesNumber(false));
    // get my petsc rank
    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);

    // test vtk output
    PETSc_Vec x, y, z, xcoor, ycoor, zcoor;  //
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, n_linear_pnts, PETSC_DECIDE);
    VecSetFromOptions(x);  //
    // get range of local variables
    VecGetOwnershipRange(x, &low, &high);
    VecGetLocalSize(x, &count);
    VecGetSize(x, &nn);
    // first duplicate x vector
    VecDuplicate(x, &y);
    VecDuplicate(x, &z);
    // get local part of vectors
    VecGetArray(x, &xp);
    VecGetArray(y, &yp);
    VecGetArray(z, &zp);

    // write coordinates

    PetscViewerASCIIPrintf(viewer, "POINTS %d double\n", nn);

    const std::vector<MeshLib::CNode*> pointVector = mesh->getNodeVector();

    for (size_t i = 0; i < (size_t)count; i++)
    {
        double const* const coords(pointVector[i]->getData());
        // now fill the vectors
        // get the pointer to current node;
        // copy local coordinates to pointer
        xp[i] = coords[0];
        yp[i] = coords[1];
        zp[i] = coords[2];
    }

    // create a sequential vector and scatter the coordinates
    // to this seq. vector on proc 0
    VecScatterCreateToZero(x, &ctx, &xcoor);
    VecScatterBegin(ctx, x, xcoor, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, x, xcoor, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterCreateToZero(y, &ctx, &ycoor);
    VecScatterBegin(ctx, y, ycoor, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, y, ycoor, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterCreateToZero(z, &ctx, &zcoor);
    VecScatterBegin(ctx, z, zcoor, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, z, zcoor, INSERT_VALUES, SCATTER_FORWARD);
    // now we have the global vector on rank 0 and can write
    if (myrank == 0)
    {
        VecGetArray(xcoor, &xp);
        VecGetArray(ycoor, &yp);
        VecGetArray(zcoor, &zp);
        for (i = 0; i < nn; i++)
        {
            PetscViewerASCIIPrintf(
                viewer, "%lg %lg %lg \n", xp[i], yp[i], zp[i]);
        }
    }

    return;
}

void LegacyVtkInterface::WriteVTKCellDataPETSC(PetscViewer viewer) const
{
    size_t nelocal = _mesh->ele_vector.size();  // local size

    // count overall length of element vector
    long numAllPoints = 0;
    for (size_t i = 0; i < nelocal; i++)
    {
        MeshLib::CElem* ele = _mesh->ele_vector[i];
        numAllPoints = numAllPoints + (ele->GetNodesNumber(false)) + 1;
    }

    PetscScalar *ep, *eglobp;  // used for pointer
    PetscInt low, high, neglob, nn;
    int i;
    VecScatter ctx;

    // get my petsc rank
    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);

    // test vtk output
    PETSc_Vec e, eglob, x;  //

    // vector for elements ...contains number of nodes connected to form a
    // element in first entry and then the node numbers
    VecCreate(PETSC_COMM_WORLD, &e);
    VecSetSizes(
        e,
        numAllPoints,
        PETSC_DECIDE);  // nummAllPoints is local lenght of the vector we need
    VecSetFromOptions(e);  //
    // get range of local variables
    VecGetOwnershipRange(e, &low, &high);
    VecGetSize(e, &neglob);
    // get local part of vector
    VecGetArray(e, &ep);

    // in order to get a global mesh the node numbers have to be corrected..
    // we have to get the position in the global node vector in order to  get
    // the correct  additive value is there a simpler way to do it (as below)?
    MeshLib::CFEMesh* mesh = fem_msh_vector[0];
    //    const int nn = mesh->getNumNodesGlobal(); //global number of nodes
    //    without shadow nodes
    const size_t n_linear_pnts(mesh->GetNodesNumber(false));

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, n_linear_pnts, PETSC_DECIDE);
    VecSetFromOptions(x);  //
    // get range of local variables
    VecGetOwnershipRange(x, &low, &high);
    VecGetSize(x, &nn);
    //  cout << " " << low << " " << high << " \n";
    // low is now the first node number

    // now we can fill the element vector
    size_t anz_elem = 0;  // here a local value
    for (size_t i = 0; i < nelocal; i++)
    {
        MeshLib::CElem* ele = _mesh->ele_vector[i];
        //	cout << ele->GetNodesNumber(false)<< " " ;

        switch (
            ele->GetElementType())  // first entry: type of element as defined
                                    // for vtk ..better than number of nodes
        {
            case MshElemType::LINE:
                ep[anz_elem] = 3;
                break;
            case MshElemType::QUAD:
                ep[anz_elem] = 9;
                break;
            case MshElemType::HEXAHEDRON:
                ep[anz_elem] = 12;
                break;
            case MshElemType::TRIANGLE:
                ep[anz_elem] = 5;
                break;
            case MshElemType::TETRAHEDRON:
                ep[anz_elem] = 10;
                break;
            case MshElemType::PRISM:  // VTK_WEDGE
                ep[anz_elem] = 13;
                break;
            case MshElemType::PYRAMID:
                ep[anz_elem] = 14;
                break;
            default:
                cerr << "PETSC VTK output::WriteVTKElementData MshElemType not "
                        "recogniced"
                     << "\n";
                break;
        }

        //   cout << ele->GetNodesNumber(false)<< " " ;
        for (size_t j = 0; j < ele->GetNodesNumber(false); j++)
        {
            ep[anz_elem + j + 1] =
                ele->GetNodeIndeces()[j] + low;  // this are the nodes
            //		cout << "DEBUG " <<  ele->GetNodeIndeces()[j] << " " ;
        }
        // cout << " \n";
        anz_elem = anz_elem + (ele->GetNodesNumber(false)) + 1;
    }

    // now scatter the vectors to process rank 0
    VecScatterCreateToZero(e, &ctx, &eglob);
    VecScatterBegin(ctx, e, eglob, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, e, eglob, INSERT_VALUES, SCATTER_FORWARD);
    // count number of elements
    if (myrank == 0)
    {
        VecGetArray(eglob, &eglobp);
        anz_elem = 0;  // here this es global number of elements
        i = 0;
        while (i < neglob)
        {
            anz_elem += 1;  // add one for the element
            int mswitch = (int)eglobp[i];
            switch (mswitch)  // first entry: type of element as defined for vtk
                              // ..better than number of nodes
            {
                case 3:  // LINE
                    i += 2 + 1;
                    break;
                case 9:  // QUAD
                    i += 4 + 1;
                    break;
                case 12:  // HEXAHEDRON:
                    i += 6 + 1;
                    break;
                case 5:  // TRIANGLE:
                    i += 3 + 1;
                    break;
                case 10:  // TETRAHEDRON:
                    i += 4 + 1;
                    break;
                case 13:  // PRISM: // VTK_WEDGE
                    i += 6 + 1;
                    break;
                case 14:  // PYRAMID:
                    i += 5 + 1;
                    break;
                default:
                    cerr << "PETSC VTK Output 1::WriteVTKElementData "
                            "MshElemType not handled, i. anz_elem. ep[i] type "
                         << i << " " << anz_elem << " " << eglobp[i] << "\n";
                    exit(1);
                    break;
            }
        }

        // write elements
        //	vtk_file << "CELLS " << numCells << " " << numAllPoints << "\n";
        PetscViewerASCIIPrintf(viewer, "CELLS %d %d\n", anz_elem, neglob);
        i = 0;
        while (i < neglob)
        {
            switch ((int)eglobp[i])  // first entry: type of element as defined
                                     // for vtk ..better than number of nodes
            {
                case 3:  // LINE
                    PetscViewerASCIIPrintf(viewer, " 2 ");
                    for (size_t j = 0; j < 2; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 2 + 1;
                    break;
                case 9:  // QUAD
                    PetscViewerASCIIPrintf(viewer, " 4 ");
                    for (size_t j = 0; j < 4; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 4 + 1;
                    break;
                case 12:  // HEXAHEDRON:
                    PetscViewerASCIIPrintf(viewer, " 6 ");
                    for (size_t j = 0; j < 6; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 6 + 1;
                    break;
                case 5:  // TRIANGLE:
                    PetscViewerASCIIPrintf(viewer, " 3 ");
                    for (size_t j = 0; j < 3; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 3 + 1;
                    break;
                case 10:  // TETRAHEDRON:
                    PetscViewerASCIIPrintf(viewer, " 4 ");
                    for (size_t j = 0; j < 4; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 4 + 1;
                    break;
                case 13:  // PRISM: // VTK_WEDGE
                    PetscViewerASCIIPrintf(viewer, " 6 ");
                    for (size_t j = 0; j < 6; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 6 + 1;
                    break;
                case 14:  // PYRAMID:
                    PetscViewerASCIIPrintf(viewer, " 5 ");
                    for (size_t j = 0; j < 5; j++)
                        PetscViewerASCIIPrintf(
                            viewer, " %d ", (long)eglobp[i + j + 1]);
                    i += 5 + 1;
                    break;
                default:
                    cerr << "PETSC VTK Output 2::WriteVTKElementData "
                            "MshElemType not handled"
                         << "\n";
                    break;
            }
            PetscViewerASCIIPrintf(viewer, " \n ");
        }

        // write cell types
        //	vtk_file << "CELL_TYPES " << numCells << "\n";
        PetscViewerASCIIPrintf(viewer, "CELL_TYPES %d \n", anz_elem);
        i = 0;
        while (i < neglob)
        {
            switch ((int)eglobp[i])  // first entry: type of element as defined
                                     // for vtk ..better than number of nodes
            {
                case 3:  // LINE
                    PetscViewerASCIIPrintf(viewer, " 3 \n");
                    i += 2 + 1;
                    break;
                case 9:  // QUAD
                    PetscViewerASCIIPrintf(viewer, " 9 \n");
                    i += 4 + 1;
                    break;
                case 12:  // HEXAHEDRON:
                    PetscViewerASCIIPrintf(viewer, " 12 \n");
                    i += 6 + 1;
                    break;
                case 5:  // TRIANGLE:
                    PetscViewerASCIIPrintf(viewer, " 5 \n");
                    i += 3 + 1;
                    break;
                case 10:  // TETRAHEDRON:
                    PetscViewerASCIIPrintf(viewer, " 10 \n");
                    i += 4 + 1;
                    break;
                case 13:  // PRISM: // VTK_WEDGE
                    PetscViewerASCIIPrintf(viewer, " 13 \n");
                    i += 6 + 1;
                    break;
                case 14:  // PYRAMID:
                    PetscViewerASCIIPrintf(viewer, " 14 \n");
                    ;
                    i += 5 + 1;
                    break;
                default:
                    cerr << "COutput::WriteVTKElementData MshElemType not "
                            "handled"
                         << "\n";
                    break;
            }
            //	PetscViewerASCIIPrintf(viewer," \n ");
        }

    }  // end myrank == 0
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2013 kg44 Implementation

   this routine is not a exact copy of the corresponding routine for !USE_PETSC
   the functionality in terms of array data output needs to be done
**************************************************************************/
void LegacyVtkInterface::WriteVTKDataArraysPETSC(PetscViewer viewer) const
{
    PetscScalar* xp;  // used for pointer
    PetscInt low, high;
    PetscInt count;

    MeshLib::CFEMesh* _mesh = fem_msh_vector[0];
    long numNodes = _mesh->GetNodesNumber(false);

    // const int nn = mesh->getNumNodesGlobal(); //global number of nodes
    // ..without shadow nodes cout << "DEBUG nNodes nn" << nNodes << " " << nn
    // <<" \n"; test vtk output
    PETSc_Vec x;  //
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, numNodes, PETSC_DECIDE);
    VecSetFromOptions(x);  //
    // get range of local variables
    VecGetOwnershipRange(x, &low, &high);
    VecGetLocalSize(x, &count);
    // get local part of vectors
    VecGetArray(x, &xp);

    // NODAL DATA
    //	vtk_file << "POINT_DATA " << numNodes << "\n";
    const size_t numPointArrays = _pointArrayNames.size();

    for (size_t k = 0; k < numPointArrays; k++)
    {
        bool toNext = false;

        // Write X, Y and Z arrays as vectors
        // KG44 needs to be re-done ...no array output so far
        {
            if (_pointArrayNames[k].find("_X") != string::npos &&
                _pointArrayNames[k + 1].find("_Y") != string::npos)
            {
                toNext = true;
            }
            // Write tensors as Eigenvectors
            // XX, XY, YY, ZZ, XZ, YZ must be present in that order
        }

        // print normal scalar fields
        if (!toNext)
        {
            string arrayName = _pointArrayNames[k];
            CRFProcess* pcs = PCSGet(arrayName, true);
            if (!pcs)
            {
            }
            else
            {
                const char* carrayName = arrayName.c_str();
                int indexDataArray = pcs->GetNodeValueIndex(arrayName);

                PetscObjectSetName((PetscObject)x, carrayName);

                for (long j = 0; j < count; j++)
                    xp[j] = pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(),
                                              indexDataArray);
                VecView(x, viewer);
            }
        }
    }

// here is the place to add the GEMS node data
#ifdef GEM_REACT
    m_vec_GEM->WriteVTKGEMValuesPETSC(
        viewer);  // kg44 export GEM internal variables like speciateion vector
                  // , phases
// etc
#endif

    // ELEMENT DATA
    // first create a PETSC vector for storing the data
    PetscScalar* ep;  // used for pointer

    //        MeshLib::CFEMesh *_mesh = fem_msh_vector[0];
    long numElem = _mesh->ele_vector.size();

    // const int nn = mesh->getNumNodesGlobal(); //global number of nodes
    // ..without shadow nodes cout << "DEBUG nNodes nn" << nNodes << " " << nn
    // <<" \n"; test vtk output
    PETSc_Vec e;  //
    VecCreate(PETSC_COMM_WORLD, &e);
    VecSetSizes(e, numElem, PETSC_DECIDE);
    VecSetFromOptions(e);  //
    // get range of local variables
    VecGetOwnershipRange(e, &low, &high);  // reuse low and high from before
    VecGetLocalSize(e, &count);            // reuse count
    // get local part of vectors
    VecGetArray(e, &ep);

    // from now on we write cell data
    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK_CELL);

    // ---------------------------------------------------------------------
    if (!_cellArrayNames.empty())
    {
        CRFProcess* pcs = this->GetPCS_ELE(_cellArrayNames[0]);

        std::vector<int> ele_value_index_vector(_cellArrayNames.size());
        if (_cellArrayNames[0].size() > 0)
            for (size_t i = 0; i < _cellArrayNames.size(); i++)
                ele_value_index_vector[i] =
                    pcs->GetElementValueIndex(_cellArrayNames[i]);

        //....................................................................
        //
        for (size_t k = 0; k < _cellArrayNames.size(); k++)
        {
            // JTARON 2010, "VELOCITY" should only write as vector, scalars
            // handled elswhere
            if (_cellArrayNames[k].compare("VELOCITY") == 0)
            {
                // kg44 to be done		vtk_file << "VECTORS velocity double "
                // <<
                // "\n"; 		this->WriteELEVelocity(vtk_file); //WW/OK
            }
            // PRINT CHANGING (OR CONSTANT) PERMEABILITY TENSOR?   // JTARON
            // 2010
            else if (_cellArrayNames[k].compare("PERMEABILITY") == 0)
            {
                /*  kg44 to be done     		vtk_file << "TENSORS
                   permeability double " << endl; for (long j = 0; j < (long)
                   _mesh->ele_vector.size(); j++)
                                        {
                                            MeshLib::CElem* ele =
                   _mesh->ele_vector[j]; CMediumProperties* MediaProp =
                                                    mmp_vector[ele->GetPatchIndex()];
                                            // KG44 22.2.2013 this is not
                   working as expected...we need to differenciate for type of
                   permeability_tensor for (size_t i = 0; i < 9; i++) vtk_file
                   << MediaProp->PermeabilityTensor(j)[i] << " ";
                        //KG44 this is buggy
                   MediaProp->PermeabilityTensor(j)[i * 3 + i] << " "; vtk_file
                   << "\n";
                                        }
                                    */
            }
            else if (ele_value_index_vector[k] > -1)
            {
                const char* carrayName = _cellArrayNames[k].c_str();
                PetscObjectSetName((PetscObject)e, carrayName);
                for (long j = 0; j < count; j++)
                    ep[j] = pcs->GetElementValue(j, ele_value_index_vector[k]);
                VecView(e, viewer);
            }
        }
        //--------------------------------------------------------------------
        ele_value_index_vector.clear();
    }
    //======================================================================
    // MAT data
    if (!_materialPropertyArrayNames.empty())
    {
        int mmp_id = -1;
        if (_materialPropertyArrayNames[0].compare("POROSITY") == 0)
            mmp_id = 0;
        // Let's say porosity
        // write header for cell data
        //            if (!wroteAnyEleData)
        //                vtk_file << "CELL_DATA " << _mesh->ele_vector.size()
        //                << "\n";
        const char* carrayName = _materialPropertyArrayNames[0].c_str();
        PetscObjectSetName((PetscObject)e, carrayName);
        for (size_t i = 0; i < (size_t)count; i++)
        {
            MeshLib::CElem* ele = _mesh->ele_vector[i];
            switch (mmp_id)
            {
                case 0:
                    ep[i] = mmp_vector[ele->GetPatchIndex()]->Porosity(i, 0.0);
                    break;
                default:
                    cout << "COutput::WriteVTKValues: no MMP values specified"
                         << "\n";
                    break;
            }
        }
        VecView(e, viewer);
    }
    // PCH: Material groups from .msh just for temparary purpose
    if (mmp_vector.size() > 1)
    {
        // write header for cell data
        //            if (!wroteAnyEleData)				//NW: check whether the
        //            header has been already written
        //                vtk_file << "CELL_DATA " << _mesh->ele_vector.size()
        //                << "\n";

        const char* carrayName = "MatGroup";
        PetscObjectSetName((PetscObject)e, carrayName);

        //            vtk_file << "SCALARS " << "MatGroup" << " int 1" << "\n";
        //            vtk_file << "LOOKUP_TABLE default" << "\n";
        for (size_t i = 0; i < (size_t)count; i++)
        {
            MeshLib::CElem* ele = _mesh->ele_vector[i];
            ep[i] = ele->GetPatchIndex();
            //                vtk_file << ele->GetPatchIndex() << "\n";
        }
        VecView(e, viewer);

        const char* carrayNameD = "Domain";
        PetscObjectSetName((PetscObject)e, carrayNameD);

        //            vtk_file << "SCALARS " << "MatGroup" << " int 1" << "\n";
        //            vtk_file << "LOOKUP_TABLE default" << "\n";
        int myrank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
        for (size_t i = 0; i < (size_t)count; i++)
        {
            // MeshLib::CElem* ele = _mesh->ele_vector[i];
            ep[i] = myrank;
            //                vtk_file << ele->GetPatchIndex() << "\n";
        }
        VecView(e, viewer);
    }

    return;
}

#else  // this is the default if petsc is not defined

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2006 KG44 Implementation
   09/2006 KG44 Output for MPI - correct OUTPUT not yet implemented
   12/2008 NW Remove ios::app, Add PCS name to VTK file name
**************************************************************************/
void LegacyVtkInterface::WriteDataVTK(int number,
                                      double simulation_time,
                                      std::string baseFilename) const
{
    baseFilename = pathJoin(defaultOutputPath, pathBasename(baseFilename));

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
    cout << "Process " << myrank << " in WriteDataVTK"
         << "\n";
#endif

    if (!_mesh)
    {
        cout << "Warning in LegacyVtkInterface::WriteVTKNodes - no MSH data"
             << "\n";
        return;
    }

    if (_processType.compare("INVALID_PROCESS") != 0)
        baseFilename += "_" + _processType;

    stringstream ss;
    // setw(4) sets the number of digits to be used
    // and creates leading zeros if necessary
    ss << setw(4) << setfill('0') << number;
    baseFilename += ss.str();
    baseFilename += ".vtk";

    // LB if(!_new_file_opened) remove(baseFilename.c_str());
    fstream vtk_file(baseFilename.c_str(), ios::out);
    vtk_file.setf(ios::scientific, ios::floatfield);
    vtk_file.precision(12);
    if (!vtk_file.good())
    {
        cout << "Could not open file for writing: " << baseFilename << "\n";
        return;
    }
    vtk_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
    // kg44 buffer the output
    char mybuffer[MY_IO_BUFSIZE * MY_IO_BUFSIZE];
    vtk_file.rdbuf()->pubsetbuf(mybuffer, MY_IO_BUFSIZE * MY_IO_BUFSIZE);
//
#endif

    this->WriteVTKHeader(vtk_file, number, simulation_time);
    this->WriteVTKPointData(vtk_file);
    this->WriteVTKCellData(vtk_file);
    this->WriteVTKDataArrays(vtk_file);
    vtk_file.close();
}

void LegacyVtkInterface::WriteVTKHeader(fstream& vtk_file,
                                        int time_step_number,
                                        double simulation_time) const
{
    vtk_file << "# vtk DataFile Version 3.0"
             << "\n";
    vtk_file << "Unstructured Grid from OpenGeoSys"
             << "\n";
    vtk_file << "ASCII"
             << "\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID"
             << "\n";

    // time information
    // see http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    vtk_file << "FIELD TimesAndCycles 2"
             << "\n";
    vtk_file << "TIME 1 1 double"
             << "\n";
    vtk_file << simulation_time << "\n";
    vtk_file << "CYLCE 1 1 long"
             << "\n";
    vtk_file << time_step_number << "\n";
}

void LegacyVtkInterface::WriteVTKPointData(fstream& vtk_file) const
{
    const std::vector<MeshLib::CNode*> pointVector = _mesh->getNodeVector();
    const size_t n_linear_pnts(_mesh->GetNodesNumber(false));
    vtk_file << "POINTS " << n_linear_pnts << " double"
             << "\n";

    for (size_t i = 0; i < n_linear_pnts; i++)
    {
        double const* const coords(pointVector[i]->getData());
        vtk_file << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
}

void LegacyVtkInterface::WriteVTKCellData(fstream& vtk_file) const
{
    size_t numCells = _mesh->ele_vector.size();

    // count overall length of element vector
    long numAllPoints = 0;
    for (size_t i = 0; i < numCells; i++)
    {
        MeshLib::CElem* ele = _mesh->ele_vector[i];
        numAllPoints = numAllPoints + (ele->GetNodesNumber(false)) + 1;
    }

    // write elements
    vtk_file << "CELLS " << numCells << " " << numAllPoints << "\n";
    for (size_t i = 0; i < numCells; i++)
    {
        MeshLib::CElem* ele = _mesh->ele_vector[i];

        // Write number of points per cell
        switch (ele->GetElementType())
        {
            case MshElemType::LINE:
                vtk_file << "2";
                break;
            case MshElemType::QUAD:
                vtk_file << "4";
                break;
            case MshElemType::HEXAHEDRON:
                vtk_file << "8";
                break;
            case MshElemType::TRIANGLE:
                vtk_file << "3";
                break;
            case MshElemType::TETRAHEDRON:
                vtk_file << "4";
                break;
            case MshElemType::PRISM:
                vtk_file << "6";
                break;
            case MshElemType::PYRAMID:
                vtk_file << "5";
                break;
            default:
                cerr << "COutput::WriteVTKElementData MshElemType not handled"
                     << "\n";
                break;
        }

        for (size_t j = 0; j < ele->GetNodesNumber(false); j++)
            vtk_file << " " << ele->getNodeIndices()[j];

        vtk_file << "\n";
    }
    vtk_file << "\n";

    // write cell types
    vtk_file << "CELL_TYPES " << numCells << "\n";

    for (size_t i = 0; i < numCells; i++)
    {
        MeshLib::CElem* ele = _mesh->ele_vector[i];

        // Write vtk cell type number (see vtkCellType.h)
        switch (ele->GetElementType())
        {
            case MshElemType::LINE:
                vtk_file << "3"
                         << "\n";
                break;
            case MshElemType::QUAD:
                vtk_file << "9"
                         << "\n";
                break;
            case MshElemType::HEXAHEDRON:
                vtk_file << "12"
                         << "\n";
                break;
            case MshElemType::TRIANGLE:
                vtk_file << "5"
                         << "\n";
                break;
            case MshElemType::TETRAHEDRON:
                vtk_file << "10"
                         << "\n";
                break;
            case MshElemType::PRISM:  // VTK_WEDGE
                vtk_file << "13"
                         << "\n";
                break;
            case MshElemType::PYRAMID:
                vtk_file << "14"
                         << "\n";
                break;
            default:
                cerr << "COutput::WriteVTKElementData MshElemType not handled"
                     << "\n";
                break;
        }
    }
    vtk_file << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2006 kg44 Implementation
   10/2006 WW Output secondary variables
   08/2008 OK MAT values
   06/2009 WW/OK WriteELEVelocity for different coordinate systems
   11/2012 WW   Rewrite this fucntion in order to have a correct vec/tensor
output
**************************************************************************/
void LegacyVtkInterface::WriteVTKDataArrays(fstream& vtk_file) const
{
    long numNodes = _mesh->GetNodesNumber(false);

    // For components of vectors and tensors. //11.2012. WW
    const int space_dim = _mesh->GetMaxElementDim();
    int vec_val_idx[3];
    int tensor_val_idx[6];

    int tensor_com = 4;

    if (space_dim == 3)
        tensor_com = 6;

    // NODAL DATA
    vtk_file << "POINT_DATA " << numNodes << "\n";
    const size_t numPointArrays = _pointArrayNames.size();

    for (size_t k = 0; k < numPointArrays; k++)
    {
        // Write tensors
        // XX, XY, YY, ZZ, XZ, YZ must be present in that order
        if (_pointArrayNames[k].find("_XX") != string::npos)
        {
            string arrayName = _pointArrayNames[k];
            CRFProcess* pcs = PCSGet(arrayName, true);
            if (!pcs)
                continue;

            tensor_val_idx[0] = pcs->GetNodeValueIndex(arrayName);
            for (int kk = 1; kk < tensor_com; kk++)
            {
                tensor_val_idx[kk] = tensor_val_idx[0] + kk;
            }

            std::string fieldName = arrayName.substr(0, arrayName.size() - 3);
            std::cout << "Writing VTK tensor: " << fieldName << "\n";
            vtk_file << "TENSORS " << fieldName << " double"
                     << "\n";

            for (long j = 0l; j < numNodes; j++)
            {
                const long node_id = _mesh->nod_vector[j]->GetIndex();
                double vector6[6];
                for (size_t component = 0;
                     component < static_cast<size_t>(tensor_com);
                     ++component)
                    vector6[component] =
                        pcs->GetNodeValue(node_id, tensor_val_idx[component]);

                vtk_file << vector6[0] << " " << vector6[1] << " " << vector6[4]
                         << "\n";
                vtk_file << vector6[1] << " " << vector6[2] << " " << vector6[5]
                         << "\n";
                vtk_file << vector6[4] << " " << vector6[5] << " " << vector6[3]
                         << "\n";
            }

            k += tensor_com - 1;
        }
        else if (_pointArrayNames[k].find("_X") != string::npos &&
                 (k + 1 < numPointArrays &&
                  _pointArrayNames[k + 1].find("_Y") != string::npos))
        {
            string arrayName = _pointArrayNames[k];
            CRFProcess* pcs = PCSGet(arrayName, true);
            if (!pcs)
                continue;

            if (pcs->getProcessType() ==
                FiniteElement::FLUID_MOMENTUM)  // 23.11.2012. WW
            {
                vec_val_idx[0] = pcs->GetNodeValueIndex(arrayName) + 1;
                for (int kk = 1; kk < space_dim; kk++)
                {
                    vec_val_idx[kk] = vec_val_idx[kk - 1] + 2;
                }
            }
            else
            {
                vec_val_idx[0] = pcs->GetNodeValueIndex(arrayName);
                for (int kk = 1; kk < space_dim; kk++)
                {
                    vec_val_idx[kk] = vec_val_idx[0] + kk;
                }
            }

            std::cout << "Writing VTK vector: " << arrayName << "\n";
            vtk_file << "VECTORS " << arrayName << " double"
                     << "\n";

            for (long j = 0l; j < numNodes; j++)
            {
                const long node_id = _mesh->nod_vector[j]->GetIndex();
                for (int kk = 0; kk < space_dim; kk++)
                {
                    vtk_file << pcs->GetNodeValue(node_id, vec_val_idx[kk])
                             << " ";
                }

                for (int kk = space_dim; kk < 3; kk++)
                {
                    vtk_file << "0.0  ";
                }

                vtk_file << "\n";
            }

            k += space_dim - 1;
        }
        else
            printScalarArray(_pointArrayNames[k], vtk_file);
    }
    //======================================================================
    // Saturation 2 for 1212 pp - scheme. 01.04.2009. WW
    // ---------------------------------------------------------------------
    CRFProcess* pcs = NULL;
    if (!_pointArrayNames.empty())  // SB added
        pcs = PCSGet(_pointArrayNames[0], true);
    if (pcs && pcs->type == 1212)
    {
        size_t i = pcs->GetNodeValueIndex("SATURATION1", true);  // JT: Latest
        vtk_file << "SCALARS SATURATION2 double 1"
                 << "\n";
        vtk_file << "LOOKUP_TABLE default"
                 << "\n";
        for (long j = 0l; j < numNodes; j++)
        {
            double val_n =
                pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(), i);
            vtk_file << 1.0 - val_n << "\n";
        }
    }
// kg44 GEM node data
#ifdef GEM_REACT
    m_vec_GEM->WriteVTKGEMValues(
        vtk_file);  // kg44 export GEM internal variables like speciateion
                    // vector , phases etc
#endif
    // ELEMENT DATA
    // ---------------------------------------------------------------------
    bool wroteAnyEleData = false;  // NW
    if (!_cellArrayNames.empty())
    {
        CRFProcess* pcs = this->GetPCS_ELE(_cellArrayNames[0]);

        std::vector<int> ele_value_index_vector(_cellArrayNames.size());
        if (_cellArrayNames[0].size() > 0)
            for (size_t i = 0; i < _cellArrayNames.size(); i++)
                ele_value_index_vector[i] =
                    pcs->GetElementValueIndex(_cellArrayNames[i]);

        vtk_file << "CELL_DATA " << (long)_mesh->ele_vector.size() << "\n";
        wroteAnyEleData = true;
        //....................................................................
        //
        for (size_t k = 0; k < _cellArrayNames.size(); k++)
        {
            // JTARON 2010, "VELOCITY" should only write as vector, scalars
            // handled elswhere
            if (_cellArrayNames[k].compare("VELOCITY") == 0)
            {
                vtk_file << "VECTORS velocity double "
                         << "\n";
                this->WriteELEVelocity(vtk_file);  // WW/OK
            }
            // PRINT CHANGING (OR CONSTANT) PERMEABILITY TENSOR?   // JTARON
            // 2010
            else if (_cellArrayNames[k].compare("PERMEABILITY") == 0)
            {
                vtk_file << "TENSORS permeability double " << endl;
                for (long j = 0; j < (long)_mesh->ele_vector.size(); j++)
                {
                    MeshLib::CElem* ele = _mesh->ele_vector[j];
                    CMediumProperties* MediaProp =
                        mmp_vector[ele->GetPatchIndex()];
                    // KG44 22.2.2013 this is not working as expected...we need
                    // to differenciate for type of permeability_tensor
                    for (size_t i = 0; i < 9; i++)
                        vtk_file << MediaProp->PermeabilityTensor(j)[i] << " ";
                    // KG44 this is buggy
                    // MediaProp->PermeabilityTensor(j)[i * 3 + i] << " ";
                    vtk_file << "\n";
                }
            }
            else if (ele_value_index_vector[k] > -1)
            {
                // NOW REMAINING SCALAR DATA  // JTARON 2010, reconfig
                vtk_file << "SCALARS " << _cellArrayNames[k] << " double 1"
                         << "\n";
                vtk_file << "LOOKUP_TABLE default"
                         << "\n";
                for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
                    vtk_file
                        << pcs->GetElementValue(i, ele_value_index_vector[k])
                        << "\n";
            }
        }
        //--------------------------------------------------------------------
        ele_value_index_vector.clear();
    }
    //======================================================================
    // MAT data
    if (!_materialPropertyArrayNames.empty())
    {
        int mmp_id = -1;
        if (_materialPropertyArrayNames[0].compare("POROSITY") == 0)
            mmp_id = 0;
        // Let's say porosity
        // write header for cell data
        if (!wroteAnyEleData)
            vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << "\n";
        wroteAnyEleData = true;
        for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
        {
            MeshLib::CElem* ele = _mesh->ele_vector[i];
            double mat_value = 0.0;
            switch (mmp_id)
            {
                case 0:
                    mat_value =
                        mmp_vector[ele->GetPatchIndex()]->Porosity(i, 0.0);
                    break;
                default:
                    cout << "COutput::WriteVTKValues: no MMP values specified"
                         << "\n";
                    break;
            }
            vtk_file << mat_value << "\n";
        }
    }
    // PCH: Material groups from .msh just for temparary purpose
    if (mmp_vector.size() > 1)
    {
        // write header for cell data
        if (!wroteAnyEleData)  // NW: check whether the header has been already
                               // written
            vtk_file << "CELL_DATA " << _mesh->ele_vector.size() << "\n";
        wroteAnyEleData = true;

        vtk_file << "SCALARS "
                 << "MatGroup"
                 << " int 1"
                 << "\n";
        vtk_file << "LOOKUP_TABLE default"
                 << "\n";
        for (size_t i = 0; i < _mesh->ele_vector.size(); i++)
        {
            MeshLib::CElem* ele = _mesh->ele_vector[i];
            vtk_file << ele->GetPatchIndex() << "\n";
        }
    }
}
#endif

/**************************************************************************
   FEMLib-Method:
   06/2009 WW/OK Implementation
**************************************************************************/
inline void LegacyVtkInterface::WriteELEVelocity(fstream& vtk_file) const
{
    int vel_ind[3] = {0, 1, 2};

    // 1D
    if (_mesh->GetCoordinateFlag() / 10 == 1)
    {
        // 0 y 0
        if (_mesh->GetCoordinateFlag() % 10 == 1)
        {
            vel_ind[0] = 1;
            vel_ind[1] = 0;
        }
        // 0 0 z
        else if (_mesh->GetCoordinateFlag() % 10 == 2)
        {
            vel_ind[0] = 2;
            vel_ind[2] = 0;
        }
    }
    // 2D
    if (_mesh->GetCoordinateFlag() / 10 == 2)
    {
        // 0 y z
        if (_mesh->GetCoordinateFlag() % 10 == 1)
        {
            vel_ind[0] = 1;
            vel_ind[1] = 2;
        }
        // x 0 z
        else if (_mesh->GetCoordinateFlag() % 10 == 2)
        {
            vel_ind[0] = 0;
            vel_ind[1] = 2;
            vel_ind[2] = 1;
        }
    }

    for (long i = 0; i < (long)_mesh->ele_vector.size(); i++)
    {
        for (int k = 0; k < 3; k++)
            vtk_file << ele_gp_value[i]->Velocity(vel_ind[k], 0) << " ";
        vtk_file << "\n";
    }
}

CRFProcess* LegacyVtkInterface::GetPCS_ELE(const string& var_name) const
{
    string pcs_var_name;
    CRFProcess* pcs = NULL;
    if (_processInfo->getProcessType() != FiniteElement::INVALID_PROCESS)
        pcs = PCSGet(_processInfo->getProcessType());
    else if (_meshTypeName.size() > 0)
        pcs = PCSGet(_meshTypeName);

    if (!pcs)
        for (size_t i = 0; i < pcs_vector.size(); i++)
        {
            pcs = pcs_vector[i];
            for (int j = 0; j < pcs->pcs_number_of_evals; j++)
            {
                pcs_var_name = pcs->pcs_eval_name[j];
                if (pcs_var_name.compare(var_name) == 0)
                    return pcs;
            }
        }
    return pcs;
}

void LegacyVtkInterface::printScalarArray(string arrayName,
                                          std::fstream& vtk_file) const
{
    CRFProcess* pcs = PCSGet(arrayName, true);
    if (!pcs)
        return;

    int indexDataArray = pcs->GetNodeValueIndex(arrayName);
    long numNodes = _mesh->GetNodesNumber(false);

    std::cout << "Writing VTK scalar: " << arrayName << "\n";
    vtk_file << "SCALARS " << arrayName << " double 1"
             << "\n";
    vtk_file << "LOOKUP_TABLE default"
             << "\n";

    for (long j = 0l; j < numNodes; j++)
        vtk_file << pcs->GetNodeValue(_mesh->nod_vector[j]->GetIndex(),
                                      indexDataArray)
                 << "\n";
}

// round very small and very large numbers in order to avoid read error in
// paraview paraview can not read number with exponents bigger/smaller than
// 300/-300
double LegacyVtkInterface::RoundDoubleVTK(double MyZahl)
{
    double rnumber;
    rnumber = MyZahl;
    if (MyZahl > 1.0e200)
        rnumber = 1.0e200;
    if (MyZahl < -1.0e200)
        rnumber = -1.0e200;
    if ((MyZahl < 1.0e-200) && (MyZahl > -1.0e-200))
        rnumber = 0.0;
    return rnumber;
}
