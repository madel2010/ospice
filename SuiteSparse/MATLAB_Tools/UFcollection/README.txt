UFcollection, Copyright 2007-2012, Timothy A. Davis, http://www.suitesparse.com

UFcollection is a MATLAB toolbox for managing the UF Sparse Matrix Collection.
If you are a MATLAB user of the collection, you would not normally need to use
this toolbox.  It contains code for creating the index for the collection (the
UF_Index.mat file in the UFget package), for creating the web pages for the
collection, and creating the Matrix Market and Rutherford/Boeing versions of
the matrices.  This code is posted here primarily so that users of the
collection can see how the matrices and their statistics were generated.

This software (UFread, specifically) also allows the user to keep a single copy
of the collection for use both inside MATLAB and outside MATLAB.  The MM/ and
RB/ versions of the collection can be read into MATLAB via UFread, even without
explicitly extracting the tar files.  They can also be read by non-MATLAB
programs.  Since the whole collection is about 8GB in size (compressed, as of
Dec 2006), this can save some space.  UFread is much slower than UFget,
however.


--------------------------------------------------------------------------------
MATLAB help for the UFcollection toolbox:
--------------------------------------------------------------------------------

  UFcollection: software for managing the UF Sparse Matrix Collection

  To create the index:

    UFindex    - create the index for the UF Sparse Matrix Collection
    UFstats    - compute matrix statistics for the UF Sparse Matrix Collection

  To create the web pages:

    UFallpages - create all web pages for the UF Sparse Matrix Collection
    UFgplot    - draw a plot of the graph of a sparse matrix
    UFint      - print an integer to a string, adding commas every 3 digits
    UFlist     - create a web page index for the UF Sparse Matrix Collection
    UFlists    - create the web pages for each matrix list (group, name, etc.)
    UFlocation - URL and top-level directory of the UF Sparse Matrix Collection
    UFpage     - create web page for a matrix in UF Sparse Matrix Collection
    UFpages    - create web page for each matrix in UF Sparse Matrix Collection
    dsxy2figxy - Transform point or position from axis to figure coords

  To create the Matrix Market and Rutherford/Boeing versions of the collection:

    UFexport     - export to Matrix Market and Rutherford/Boeing formats
    UFread       - read a Problem in Matrix Market or Rutherford/Boeing format
    UFwrite      - write a Problem in Matrix Market or Rutherford/Boeing format
    UFfull_read  - read a full matrix using a subset of Matrix Market format
    UFfull_write - write a full matrix using a subset of Matrix Market format

  Example:
    UFindex       % create index (UF_Index.mat) for use by UFget
    UFallpages    % create all web pages for the UF Sparse Matrix Collection

  Requires UFget, CSparse, CHOLMOD, AMD, COLAMD, RBio, and METIS.


--------------------------------------------------------------------------------
Files:
--------------------------------------------------------------------------------

    Contents.m		    MATLAB help
    dsxy2figxy.m	    convert XY points for plot annotations
    Makefile		    Unix/Linux installation, or use UFcollection_install
    README.txt		    this file
    UFallpages.m	    create all web pages
    UFexport.m		    export to MM and RB
    UFcollection_install.m  installation
    UFfull_read.m	    read a full matrix
    UFfull_write.c	    write a full matrix
    UFfull_write.m	    MATLAB help for UFfull_write
    UFgplot.m		    plot a graph
    UFindex.m		    create UF_Index.mat
    UFint.m		    print an integer
    UFlist.m		    create a web page index
    UFlists.m		    create all web page indices
    UFlocation.m	    URL and directory for the collection
    UFpage.m		    create a web page for a matrix
    UFpages.m		    create web pages for all matrices
    UFread.m		    read a Problem
    UFstats.m		    compute statistics about a matrix
    UFwrite.m		    write a Problem

./Doc:
    gpl.txt		    GNU GPL license
    License.txt

--------------------------------------------------------------------------------
To add a matrix to the collection:
--------------------------------------------------------------------------------

These instructions are for the maintainer of the collection (that is, just
notes to myself), but they also indicate how the above software is used.

Requires most of SuiteSparse (UFget, CHOLMOD, AMD, COLAMD, CSparse, RBio, and
UFcollection), and METIS 4.0.1.

1) Get the matrix into MATLAB (method depending on how the matrix was
    submitted).  Use load and sparse2, RBread, mread, or specialized code
    written just for that matrix.

2) Add the matrix to the end of UF_Listing.txt (a line in the form Group/Name).

3) Create a new directory /cise/research/sparse/public_html/mat/Group,
    where Group is the new matrix group.  Add a README.txt file to this
    directory, the first line of which is a one-line summary that will appear
    in the top-level web page for the collection.  Skip this step if adding a
    matrix to an existing group.

4) Create the Problem struct; type "help UFwrite" for details.  Required fields:

    Problem.name    full name of the matrix (Group/Name)
    Problem.title   short descriptive title
    Problem.A	    the sparse matrix
    Problem.id	    integer corresponding to the line number in UF_Listing.txt
    Problem.date    date the matrix was created, or added to the collection
    Problem.author  matrix author
    Problem.ed	    matrix editor/collector
    Problem.kind    a string.  For a description, see below.

    optional fields:

    Problem.Zeros   binary pattern of explicit zero entries
    Problem.b	    right-hand-side
    Problem.x	    solution
    Problem.notes   a char array
    Problem.aux	    auxiliary matrices (contents are problem dependent)

    Save to a MATLAB mat-file.  In the mat directory, do:

    save (Problem.name, 'Problem', '-v7') ;

    or for very large problems (ids 1903 and 1905, for example):

    save (Problem.name, 'Problem', '-v7.3') ;

5) Compute matrix statistics and extend the UF_Index:

    UFindex (ids)

    where ids is a list of the new matrix id's.  Updated UF_Index.mat and
    UFstats.csv files are placed in the current working directory.
    Copy the new UF_Index.mat file to /cise/research/sparse/public_html/mat.
    Copy the new UFstats.csv file to /cise/research/sparse/public_html/matrices.

6) Update the web pages:

    In the /cise/research/sparse/public_html directory, do:

    UFlists
    UFpages (1, ids)

7) Export the matrix in Matrix Market and Rutherford/Boeing formats.

    UFexport (ids)

    or

    UFexport (ids, 'check')

    then tar and compress the resulting MM/Group/Name and RB/Group/Name
    directories, one per Problem (if UFexport has not already done so).
    Copy the MM and RB matrices from the 2sparse/MM and /RB directories into
    the sparse/public_html directory.

8) Make the collection world-readable.  In /cise/research/sparse/public_html do:

    chmod -R og+rX mat matrices MM RB

9) Optional: if a new group was added, manually edit the
    /cise/research/sparse/public_html/matrices/index.html file, adding a
    new thumbnail image to the Sample Gallery.  If a new Problem.kind was
    introduced, describe it in the matrices/kind.html file.


--------------------------------------------------------------------------------
Problem.kind
--------------------------------------------------------------------------------

Problems with 2D/3D geometry

    2D/3D problem
    acoustics problem
    computational fluid dynamics problem
    computer graphics/vision problem
    electromagnetics problem
    materials problem
    model reduction problem
    robotics problem
    semiconductor device problem
    structural problem
    thermal problem 

Problems that normally do not have 2D/3D geometry

    chemical process simulation problem
    circuit simulation problem
    counter-example problem: Some of these may have 2D/3D geometry.
    economic problem
    frequency-domain circuit simulation problem
    least squares problem
    linear programming problem
    optimization problem
    power network problem
    statistical/mathematical problem
    theoretical/quantum chemistry problem
    combinatorial problem

Graph problems

    This problem includes the graph or multigraph keyword. It is a network or
    graph. A graph may or may not have 2D/3D geometry (typically it does not).
    Several secondary phrases can be included in Problem.kind:

        directed or undirected: A graph is either directed, undirected, or
        bipartite. (bipartite graphs are always undirected). If not bipartite,
        the matrix will always be square. Unsymmetric permutations of the
        matrix have no meaning. If directed, the edge (i,j) is not the same as
        (j,i), and the matrix will normally be unsymmetric. If undirected, the
        edges (i,j) and (j,i) are the same, and the matrix is always symmetric.

        weighted: If the graph has edge weights, this word will appear. The
        edge weight of edge (i,j) is the value of A(i,j). This phrase is used
        for a graph only, never for a multigraph. If the graph is not weighted,
        the matrix is binary.

        bipartite: If the rows and columns of the matrix reflect different sets
        of nodes. The matrix A is normally rectangular, but can be square. Any
        permutation (unsymmetric or symmetric) of the matrix is meaningful.

        random: This Problem has been randomly generated. It is included in the
        UF Sparse Matrix Collection only because it has been used as a standard
        benchmark. Randomly generated problems are otherwise excluded from the
        collection.  multigraph or graph: If the matrix represents a
        multigraph, then A(i,j) reflects the number of edges (i,j). The edges
        themselves are always unweighted. If the matrix represents a graph,
        then A(i,j) is either 0 or 1 for an unweighted graph, or the weight of
        edge (i,j) otherwise. 

