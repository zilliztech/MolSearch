**This repo will no longer be maintained, please visit https://github.com/milvus-io/bootcamp/tree/master/solutions/molecular_similarity_search.**
# MolSearch


[English](README.md) | [中文版](CN_README.md)

# Requirements 

- ### [Milvus CPU](https://milvus.io/docs/v0.10.4/milvus_docker-cpu.md)

- ### postgres

- ### RDKit

# Deploy

## 1. Run Milvus Docker

This demo uses Milvus-0.10.4-CPU，please refer to https://milvus.io/docs/v0.10.4/milvus_docker-cpu.md。

```
# Start Milvus
$ docker run -d --name milvus_cpu_0.10.4 \
-p 19530:19530 \
-p 19121:19121 \
-v /home/$USER/milvus/db:/var/lib/milvus/db \
-v /home/$USER/milvus/conf:/var/lib/milvus/conf \
-v /home/$USER/milvus/logs:/var/lib/milvus/logs \
-v /home/$USER/milvus/wal:/var/lib/milvus/wal \
milvusdb/milvus:0.10.4-cpu-d120220-e72454
```

## 2. Run molsearch-webserver docker

```
$ docker run -td -p 35001:5000 -e "MILVUS_HOST=192.168.1.85" -e "MILVUS_PORT=19530" -e "PG_HOST=192.168.1.85" -e "PG_PORT=5432" zilliz/molsearch-webserver:0.2.0
```

The description for parameters:

| Parameter                     | Description                                                  |
| ----------------------------- | ------------------------------------------------------------ |
| -p 35001:5000                 | -p specifies pot mapping between the host and the image.     |
| -e "MILVUS_HOST=192.168.1.85" | -e specifies the system parameter mapping between the host and the image. Pease update `192.168.1.25` to the IP address of the Milvus docker. |
| -e "MILVUS_PORT=19530"        | Update `19530` to the port of Milvus docker.                 |

## 3. Run molsearch-webclient docker

```
$ docker run -td -p 8001:80 -e API_URL=http://192.168.1.85:35001  zilliz/molsearch-webclient:0.1.0
```

> Note: Please update `192.168.1.85` to the IP address of the Milvus docker.

## 4. Import data to Milvus

Import '.smi' data into Milvus, where the first column is smiles and the second column is the id number, for example:

o1c(C(O)CNC(C)(C)C)cc2c1c(CC(=O)OC(C)(C)C)ccc2    10001

```bash
# Run the following command under the script directory
$ cd script
$ python insert_data.py -f test_1w.smi
```

> The data is from  [pubchem](ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF) and [zinc](http://zinc.docking.org/tranches/home/), which extracted 10,000 pieces of data.

You can enter the IP and port where molsearch-webclient is launched in your browser.

```bash
192.168.1.85:8001
```

# How to use

MolSearch is an open source molecular search software based on [Milvus](https://github.com/milvus-io/milvus) and [MolView](https://github.com/molview/legacy), which has six main feature: editor, chemical formula, molecular search, tool classes, 3D model, Jmol tools.

![img](./pic/molsearch.png)

## Drawing structural formulas

### Top toolbar

![img](./pic/draw1.png)            

- **Trash:** clear the entire canvas
- **Eraser:** erase atoms, bonds or the current selection
- **Undo/redo:** undo or redo your recent changes
- Selection tools: all these tool can be used to drag the current selection or individual atoms and bonds. You can add/remove atoms and bonds to the selection by clicking them. If you have selected a separate fragment, you can rotate it by dragging an atom in the selection. You can delete the selection using the **DEL** key or using the eraser tool. Each tool has different behavior for the right mouse button:
  - **Drag:** move the entire molecule (you can already use the left mouse button for this)
  - **Rectangle select:** select atoms and bonds using a rectangular selection area
  - **Lasso select:** select atoms and bonds by drawing a freehand selection area
- **Color mode:** display atoms and bonds using colors
- **skeleton mode:** displays all C and H atoms instead of skeletal display
- **Center:** centers the whole molecule
- **Clean:** cleans the structural formula using an external service
- **to 3D:** converts the structural formula into a 3D model

#### Left toolbar

![](./pic/draw2.png)             

- **Bonds:** pick one of the bond types (single, double, triple, up, down) and add or modify bonds
- **Fragments:** pick one of the fragments (benzene, cyclopropane, etc.) and add fragments
- **Chain:** create a chain of carbon atoms
- **Charge:** increment (+) or decrement (-) the charge of atoms

#### Right toolbar

![](./pic/draw3.png)

In this toolbar you can select from a number of elements, you can also pick an element from the periodic table using the last button. You can use the element to create new atoms or modify existing atoms.

## Finding structures

 ![](./pic/load.png)

You can load molecules, just type what you are looking for and a list of available molecules will appear.

## Search

These functions allow you to perform some advanced searches through the database using the structural formula from the sketcher.

1. **Similarity search:** search for compounds with a similar structural formula
2. **Substructure search:** search for compounds with the current structure as subset
3. **Superstructure search:** search for compounds with the current structure as superset

## Tools

The **Tools** menu contains several utility functions which are listed below.

- **Structural formula image:** sketcher snapshot (PNG with alpha channel)
- **3D model image:** model snapshot (PNG)
- **MOL file:** exports a MDL Molfile from the 3D model **(common molecules)**

#### Information card

This collects and displays information about the structural formula.

## 3D Model

The **Model** menu contains some general functions for the 3D model.

#### Reset

This function sets the model position, zoom and rotation back to default.

#### Representation

You can choose from a list of different molecule representations including; ball and stick, stick, van der Waals spheres, wireframe and lines. Macromolecules are automatically drawn using ribbons.

#### Background

You can switch between a black, gray or white background. The default background is black (exported images from GLmol or ChemDoodle have a transparent background)

#### Engines

You can choose from three different render engines: **GLmol**, **Jmol** and **ChemDoodle**. GLmol is used as default render engine. GLmol and ChemDoodle are based on WebGL, a browser technology to support 3D graphics. If WebGL is not available in your browser, Jmol will be used for all rendering.

MolSearch automatically switches to:

1. **Jmol** if you execute functions from the Jmol menu
2. **GLmol** if you load macromolecules (due to significant higher performance)
3. **ChemDoodle** if you load a crystal structure (GLmol cannot render crystal structures)

## Jmol tools

The **Jmol** menu offers some awesome Jmol-only functions and calculations.

#### Clear

Clears all executed calculations and measurements.

#### High Quality

Enables High Quality rendering in Jmol (enabled by default on fast devices) When turned off, anti-aliasing is disabled and the model is drawn using lines while transforming it.

#### Calculations

You can perform the following Jmol calculations in Jmol:

- **MEP surface lucent/opaque:** calculates and projects molecular analysis software electrostatic potential on a translucent or opaque van der Waals surface
- **Charge:** calculates and projects atomic charge as text label and white to atom color gradient
- **Bond dipoles:** calculates and draws individual bond dipoles
- **Overall dipole:** calculates and draws net bond dipole
- **Energy minimization:** executes an interactive MMFF94 energy minimization *(note that this function only executes a maximum of 100 minimization steps at a time)*

#### Measurement

You can measure distance, angle and torsion using Jmol. You can activate and deactivate one of these measurement types via the Jmol menu.

- **Distance** distance between two atoms in nm
- **Angle** angle between two bonds in degrees
- **Torsion** torsion between four atoms in degrees

Note that in some cases, the resolved 3D model is only an approach of the real molecule, this means you have to execute an **Energy minimization** in order to do reliable measurements.
