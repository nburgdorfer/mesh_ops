/* Eric Joyce, Stevens Institute of Technology, 2019; ejoyce@stevens.edu

   19dec19:  Receives a mesh and returns a reduced mesh that has dropped all triangles with areas >= threshold.
             Assumes vertices of the given mesh have attributes x, y, z, r, g, b, and isCam.
*/

#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
#define __DROPTHRESHOLD_DEBUG 1
*/
/*
#define __DROPTHRESHOLD_LOAD_DEBUG 1
*/

typedef struct VertexType
  {
    float x;                                                        //  The specifications of these data types
    float y;                                                        //  come from the ETH PLY file headers.
    float z;                                                        //  Even if we do not use all of these for
    unsigned char r;                                                //  our purposes, we need to know how large
    unsigned char g;                                                //  they are so that we can skip past them
    unsigned char b;                                                //  in the binary data.
    bool cam;                                                       //  uchar in PLY spec: convert to bool
    bool keep;                                                      //  Marked as a Vertex to be kept after triangles have been dropped
  } Vertex;

typedef struct FaceType
  {
                                                                    //  An initial value, n, is always going to be = 3. Ignore it.
    unsigned int v[3];                                              //  Indices of Vertices making up this Face
    double barycenter[3];                                           //  X, Y, and Z components of the Face barycenter
    double normal[3];                                               //  X, Y, and Z components of unit normal to the Face
    double area;                                                    //  Area of this face
  } Face;

/*
g++ -Wall areaTh.cpp -lm -o areaTh
./areaTh delivery_area.mincut.KB.ply 1.0
*/

using namespace std;

void writePLY(const char*, Vertex**, Face*, unsigned int*, unsigned int*, double);
void writeAuxiliaryFile(char*, unsigned int*, unsigned int*);
bool loadMesh(const char*, unsigned int*, Vertex**, Face**);
bool getMeshInfo(const char*, unsigned int**);
void computeLongestLeg(Vertex*, Face**, unsigned int);
void computeFaceBarycenter(Vertex*, Face**, unsigned int);
void computeFaceArea(Vertex*, Face**, unsigned int);
void computeFaceNormal(Vertex*, Face**, unsigned int);
float crossProduct(float*, float*, float**);
void quicksort(double**, unsigned int**, unsigned int, unsigned int);
unsigned int partition(bool, double**, unsigned int**, unsigned int, unsigned int);
void usage(void);

int main(int argc, char* argv[])
  {
    Vertex* V;                                                      //  Array of vertices
    Face* F;                                                        //  Array of faces

    double threshold = 1.00;                                        //  Area threshold, above which we drop
    unsigned int* fileInfo;                                         //  Array of data about the given PLY file
    unsigned char fileInfoLen;                                      //  Length of that array

    double* areas;                                                  //  Array of Face areas sorted from SMALLEST to LARGEST
    unsigned int* indices;                                          //  Accompanying array of Face indices
    char auxFilename[128];
    ifstream fh;
    unsigned int i, j, tmp;
    bool found = false;

    fileInfoLen = 4;                                                //  There are 4 things we need to know about the PLY:
    if(argc < 2)                                                    //  [0] Number of vertices
      {                                                             //  [1] Number of faces
        usage();                                                    //  [2] Number of tetrahedra (not needed here)
        return 0;                                                   //  [3] Offset of the first byte of PLY binary data
      }                                                             //  And make sure we have enough arguments!
    if(argc > 2)                                                    //  This argument is optional:
      threshold = (double)atof(argv[2]);                            //  If given a 'threshold', use that

    i = strlen(argv[1]) - 1;                                        //  Back up until we've left the file extension
    while(argv[1][i] != '.')
      i--;
    for(j = 0; j < i; j++)
      auxFilename[j] = argv[1][j];                                  //  Copy the source file name up to the extension
    auxFilename[j] = '.';  j++;
    auxFilename[j] = 'a';  j++;
    auxFilename[j] = 'r';  j++;
    auxFilename[j] = 'e';  j++;
    auxFilename[j] = 'a';  j++;
    auxFilename[j] = '\0';                                          //  NULL-terminate
                                                                    //  Allocate file info array
    if((fileInfo = (unsigned int*)malloc(fileInfoLen * sizeof(int))) == NULL)
      {
        cout << "ERROR: Unable to allocate file info array." << endl;
        return 0;
      }

    if(!getMeshInfo(argv[1], &fileInfo))                            //  Get PLY header information
      {
        cout << "ERROR: Unable to allocate file info array." << endl;
        return 0;
      }

    if(!loadMesh(argv[1], fileInfo, &V, &F))                        //  Only load vertices and faces
      {
        if(fileInfo[0] > 0)
          free(V);
        if(fileInfo[1] > 0)
          free(F);
        free(fileInfo);
        return 0;
      }
    if((indices = (unsigned int*)malloc(fileInfo[1] * sizeof(int))) == NULL)
      {
        cout << "ERROR: Unable to allocate Face index array." << endl;
        if(fileInfo[0] > 0)
          free(V);
        if(fileInfo[1] > 0)
          free(F);
        free(fileInfo);
        free(areas);
        return 0;
      }

    fh.open(auxFilename);                                           //  Check to see whether an auxiliary file for this PLY exists!
    if(fh.good())                                                   //  Save SOOOO much time!
      {
        #ifdef __DROPTHRESHOLD_DEBUG
        cout << "FOUND AREA FILE!" << endl;
        #endif

        found = true;
        i = 0;                                                      //  Read indices in order
        while(!fh.eof() && i < fileInfo[1])
          {
            fh.read(reinterpret_cast<char*>(&tmp), sizeof(int));
            indices[i] = tmp;
            i++;
          }
        fh.close();
      }
    else
      {
                                                                    //  Allocate Face area array
        if((areas = (double*)malloc(fileInfo[1] * sizeof(double))) == NULL)
          {
            cout << "ERROR: Unable to allocate Face area array." << endl;
            if(fileInfo[0] > 0)
              free(V);
            if(fileInfo[1] > 0)
              free(F);
            free(fileInfo);
            return 0;
          }
                                                                    //  Allocate Face index array
        for(i = 0; i < fileInfo[1]; i++)                            //  Fill initial values of Face arrays
          {
            areas[i] = F[i].area;
            indices[i] = i;
          }

        quicksort(&areas, &indices, 0, fileInfo[1] - 1);            //  Sort the arrays according to area
        writeAuxiliaryFile(auxFilename, indices, fileInfo);         //  Save this information to auxiliary file!
      }

    writePLY("areaTh.ply", &V, F, fileInfo, indices, threshold);    //  Create the new file

    if(fileInfo[0] > 0)                                             //  Clean up, go home
      free(V);
    if(fileInfo[1] > 0)
      {
        free(F);
        if(!found)
          free(areas);
        free(indices);
      }
    free(fileInfo);

    return 0;
  }

/* Write the new, decimated mesh to file.
   Along the way, identify which vertices are necessary to the decimated mesh.
   Note that here we DO include the isCam attribute because this mesh is likely to go onto PSAMPLE. */
void writePLY(const char* filename, Vertex** V, Face* F, unsigned int* info, unsigned int* faceIndices, double threshold)
  {
    ofstream fh;                                                    //  Output file handle
    unsigned int cutoff = 0;                                        //  Length (not the index) of the diminished Face array

    unsigned int* map;                                              //  Lookup table: old index of Vertex ---> new index of Vertex
    unsigned int curr = 0;                                          //  Current new vertex index

    float x, y, z;
    unsigned char r, g, b, holder;

    unsigned int i, j;

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "writePLY()" << endl;
    #endif
                                                                    //  Position cutoff at the first triangle that's too large
    while(cutoff < info[1] && F[ faceIndices[cutoff] ].area < threshold)
      cutoff++;

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "  Threshold = " << F[ faceIndices[cutoff] ].area << endl;
    #endif
                                                                    //  Go through the Faces at these indices
    for(i = 0; i < cutoff; i++)                                     //  and mark all their Vertices as 'to-be-kept'.
      {
        for(j = 0; j < 3; j++)                                      //  For every Vertex of Face [ faceIndices[i] ]...
          (*V)[ F[ faceIndices[i] ].v[j] ].keep = true;             //  mark it as necessary to keep
                                                                    //  (remember, they'd all been marked false when loaded)
                                                                    //  Increment number of marked Vertices
      }                                                             //  Two separate loops takes less time than two nested loops!!!
    if((map = (unsigned int*)malloc(info[0] * sizeof(int))) == NULL)//  Allocate Vertex lookup array
      {
        cout << "ERROR: Unable to allocate vertex map array." << endl;
        exit(1);
      }
    for(i = 0; i < info[0]; i++)                                    //  Blank out lookup array with INT_MAX
      map[i] = INT_MAX;
    for(i = 0; i < info[0]; i++)                                    //  Count up all marked Vertices
      {
        if((*V)[i].keep)                                            //  If this Vertex is alive
          {
            map[i] = curr;
            curr++;
          }
      }

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "  " << +curr << " / " << +info[0] << " vertices in reduced mesh." << endl;
    cout << "  " << +cutoff << " / " << +info[1] << " faces in reduced mesh." << endl;
    #endif

    fh.open(filename);                                              //  Open the file
    fh << "ply" << endl;                                            //  Ply format header
    fh << "format binary_little_endian 1.0" << endl;                //  Declare endianness
    fh << "element vertex " << +curr << endl;                       //  Declare number of vertices
    fh << "property float x" << endl;                               //  Vertex properties are the same as in the original PLY
    fh << "property float y" << endl;
    fh << "property float z" << endl;
    fh << "property uchar red" << endl;
    fh << "property uchar green" << endl;
    fh << "property uchar blue" << endl;
    fh << "property uchar isCam" << endl;
    fh << "element face " << +cutoff << endl;                       //  Declare number of faces
    fh << "property list uchar int vertex_index" << endl;
    fh << "end_header" << endl;                                     //  Close the header

    for(i = 0; i < info[0]; i++)                                    //  Write all (surviving) Vertices
      {
        if(map[i] != INT_MAX)
          {
            x = (*V)[ i ].x;
            y = (*V)[ i ].y;
            z = (*V)[ i ].z;
            r = (*V)[ i ].r;
            g = (*V)[ i ].g;
            b = (*V)[ i ].b;
            holder = ((*V)[ i ].cam) ? 1 : 0;

            fh.write(reinterpret_cast<char*>(&x), sizeof(float));
            fh.write(reinterpret_cast<char*>(&y), sizeof(float));
            fh.write(reinterpret_cast<char*>(&z), sizeof(float));
            fh.write(reinterpret_cast<char*>(&r), sizeof(char));
            fh.write(reinterpret_cast<char*>(&g), sizeof(char));
            fh.write(reinterpret_cast<char*>(&b), sizeof(char));
            fh.write(reinterpret_cast<char*>(&holder), sizeof(char));
          }
      }

    for(i = 0; i < cutoff; i++)                                     //  Write each Face
      {
        holder = 3;
        fh.write(reinterpret_cast<char*>(&holder), sizeof(char));   //  Easy: just write 3

        j = map[ F[ faceIndices[i] ].v[0] ];
        fh.write(reinterpret_cast<char*>(&j), sizeof(int));         //  Write new index of Vertex-0

        j = map[ F[ faceIndices[i] ].v[1] ];
        fh.write(reinterpret_cast<char*>(&j), sizeof(int));         //  Write new index of Vertex-1

        j = map[ F[ faceIndices[i] ].v[2] ];
        fh.write(reinterpret_cast<char*>(&j), sizeof(int));         //  Write new index of Vertex-2
      }
    if(curr > 0)                                                    //  Release array
      free(map);

    fh.close();                                                     //  Close the file. Done!

    return;
  }

/* Write Face indices in INCREASING ORDER of area to an auxiliary file.
   That way, we only need to compute areas ONCE!
   Save SO MUCH TIME! */
void writeAuxiliaryFile(char* filename, unsigned int* indices, unsigned int* info)
  {
    ofstream fh;                                                    //  Output file handle
    unsigned int i, index;

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "writeAuxiliaryFile()" << endl;
    #endif

    fh.open((const char*)filename);                                 //  Open the file
    for(i = 0; i < info[1]; i++)
      {
        index = indices[i];
        fh.write(reinterpret_cast<char*>(&index), sizeof(int));
      }
    fh.close();                                                     //  Close the file. Done!

    return;
  }

/* Read the binary data of the given PLY file and load the given arrays.
   Return true on success.
   Return false on failure. */
bool loadMesh(const char* filename, unsigned int* info, Vertex** V, Face** F)
  {
    ifstream fh;                                                    //  File handle

    float x, y, z;                                                  //  Elements of Vertex
    unsigned char r, g, b;
    unsigned char isCam;
    unsigned char num;                                              //  Elements of Face
    unsigned int v0, v1, v2;
    unsigned char phase;                                            //  Into which array are we reading
    unsigned int i;
    bool breakout;

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "load(" << filename;
    cout << ", [" << +info[0] << ", " << +info[1] << ", " << +info[2] << ", " << +info[3] << "])" << endl;
    #endif

    fh.open(filename, ios::binary);                                 //  Open the given file (for binary reading this time)
    fh.seekg(info[3], ios::beg);                                    //  Jump past the header to the PLY data
                                                                    //  Allocate arrays for vertices, faces, tetrahedra
    if(((*V) = (Vertex*)malloc(info[0] * sizeof(Vertex))) == NULL)
      {
        cout << "ERROR: Unable to allocate vertex buffer." << endl;
        return false;
      }
    if(((*F) = (Face*)malloc(info[1] * sizeof(Face))) == NULL)
      {
        cout << "ERROR: Unable to allocate face buffer." << endl;
        free((*V));
        return false;
      }

    i = 0;
    phase = 0;
    breakout = false;
    while(!fh.eof() && !breakout)
      {
        switch(phase)
          {
            case 0:                                                 //  Read elements from file
                     fh.read(reinterpret_cast<char*>(&x), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&y), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&z), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&r), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&g), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&b), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&isCam), sizeof(char));
                     (*V)[i].x = x;                                 //  Write elements to struct in array
                     (*V)[i].y = y;
                     (*V)[i].z = z;
                     (*V)[i].r = r;
                     (*V)[i].g = g;
                     (*V)[i].b = b;
                     (*V)[i].cam = (isCam != 0) ? true : false;
                     (*V)[i].keep = false;                          //  Initialize everything to false; mark true if necessary in writePLY()
                     #ifdef __DROPTHRESHOLD_LOAD_DEBUG
                     cout << "Vertex[" << +i << "]:" << endl;
                     cout << "    x: " << (*V)[i].x << " <-- " << x << endl;
                     cout << "    y: " << (*V)[i].y << " <-- " << y << endl;
                     cout << "    z: " << (*V)[i].z << " <-- " << z << endl;
                     cout << "    r: " << +(*V)[i].r << " <-- " << +r << endl;
                     cout << "    g: " << +(*V)[i].g << " <-- " << +g << endl;
                     cout << "    b: " << +(*V)[i].b << " <-- " << +b << endl;
                     if((*V)[i].cam)
                       cout << "  cam: 1";
                     else
                       cout << "  cam: 0";
                     cout << " <-- " << +isCam << endl;
                     #endif
                     if(++i == info[0])                             //  If we've read sufficient,
                       {                                            //  then change phases and reset the index
                         phase++;
                         i = 0;
                       }
                     break;
            case 1:                                                 //  Read elements from file
                                                                    //  Read 'num' just to advance the file pointer; we don't care.
                     fh.read(reinterpret_cast<char*>(&num), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&v0), sizeof(int));
                     fh.read(reinterpret_cast<char*>(&v1), sizeof(int));
                     fh.read(reinterpret_cast<char*>(&v2), sizeof(int));
                     (*F)[i].v[0] = v0;
                     (*F)[i].v[1] = v1;
                     (*F)[i].v[2] = v2;
                     computeFaceBarycenter((*V), F, i);             //  Compute and store the barycenter for this Face
                     computeFaceNormal((*V), F, i);                 //  Compute and store the unit normal for this face
                     computeFaceArea((*V), F, i);                   //  Compute and store the area of this Face

                     #ifdef __DROPTHRESHOLD_LOAD_DEBUG
                     cout << "Face:" << endl;
                     cout << "  v0:     " << +(*F)[i].v[0] << " <-- " << +v0 << endl;
                     cout << "  v1:     " << +(*F)[i].v[1] << " <-- " << +v1 << endl;
                     cout << "  v2:     " << +(*F)[i].v[2] << " <-- " << +v2 << endl;
                     cout << "  normal: " << (*F)[i].normal[0] << ", " << (*F)[i].normal[1] << ", " << (*F)[i].normal[2] << endl;
                     cout << "  area:   " << (*F)[i].area << endl;
                     #endif
                     if(++i == info[1])                             //  If we've read sufficient,
                       breakout = true;                             //  then we're done!
                     break;
          }
      }
    fh.close();                                                     //  Close the file

    return true;
  }

/* Read the given PLY file and find out all the information we need to know about its binary data.
   Write the following details into the given array, 'info'.
     [0] Number of vertices                                         e.g.  element vertex 31480
     [1] Number of faces                                            e.g.  element face 408330
     [2] Number of tetrahedra                                       e.g.  element tetrahedra 204129
     [3] Offset of the first byte of PLY binary data
   Return true on success.
   Return false on failure. */
bool getMeshInfo(const char* filename, unsigned int** info)
  {
    ifstream fh;                                                    //  File handle
    unsigned int contentsLen;                                       //  Length of given file in bytes
    char* buffer;                                                   //  Contains file contents
    char line[256];                                                 //  Line buffer
    unsigned int i, j;
    unsigned int linectr;
    unsigned char spctr;                                            //  Count spaces
    char* word0;                                                    //  Substring: first word found in line
    char* word1;                                                    //  Substring: second word found in line
    char* word2;                                                    //  Substring: third word found in line
    unsigned char len0, len1, len2;                                 //  Lengths of these words

    #ifdef __DROPTHRESHOLD_DEBUG
    cout << "getMeshInfo(" << filename << ")" << endl;
    #endif

    fh.open(filename);                                              //  Open the given file
    fh.seekg(0, ios::end);                                          //  Jump to the end
    contentsLen = fh.tellg();                                       //  Find out how many bytes we need to allocate
    fh.seekg(0, ios::beg);                                          //  Jump back to the head of the file
    if((buffer = (char*)malloc(contentsLen * sizeof(char))) == NULL)//  Allocate the bytes
      return 0;
    fh.read(buffer, contentsLen);                                   //  Read in entire file
    fh.close();                                                     //  Close the file

    i = 0;                                                          //  We're only interested in the ASCII header for now,
    linectr = 0;                                                    //  so read line by line until "end_header" is found.
    while(i < contentsLen)                                          //  (This is an indulgent bound, but all we have so far.)
      {
        if(buffer[i] == '\n')                                       //  End of a line: time to examine the line buffer
          {
            line[linectr] = '\0';                                   //  Null-terminate the line buffer

            spctr = 0;                                              //  How many spaces on this header line?
            for(j = 0; j < linectr; j++)                            //  Count spaces
              {
                if(line[j] == ' ')
                  spctr++;
              }
            if(spctr == 2)                                          //  Two spaces found means three words on this line: a candidate!
              {
                j = 0;
                spctr = 0;                                          //  Reset space counter
                len0 = 0;                                           //  Set all initial lengths to zero
                len1 = 0;
                len2 = 0;
                while(j <= linectr)
                  {
                    switch(spctr)                                   //  Every iteration adds a character
                      {                                             //  to one of our three words.
                        case 0:  if(++len0 == 1)                    //  Expand the array
                                   {
                                     if((word0 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word0 = (char*)realloc(word0, len0 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                                 break;
                        case 1:  if(++len1 == 1)                    //  Expand the array
                                   {
                                     if((word1 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word1 = (char*)realloc(word1, len1 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                                 break;
                        default: if(++len2 == 1)                    //  Expand the array
                                   {
                                     if((word2 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word2 = (char*)realloc(word2, len2 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                      }

                    if(line[j] == ' ')                              //  Found a space: one word ends; point to the next.
                      {
                        switch(spctr)
                          {
                            case 0:  word0[len0 - 1] = '\0';  break;//  Cap word-0
                            case 1:  word1[len1 - 1] = '\0';  break;//  Cap word-1
                          }
                        spctr++;                                    //  Increment space counter
                      }
                    else if(line[j] == '\0')                        //  Found a null
                      word2[len2 - 1] = '\0';                       //  Cap word-2
                    else                                            //  Found a character
                      {                                             //  The current array has already been expanded,
                        switch(spctr)                               //  So just add the current character to the current array!
                          {
                            case 0:  word0[len0 - 1] = line[j];  break;
                            case 1:  word1[len1 - 1] = line[j];  break;
                            default: word2[len2 - 1] = line[j];  break;
                          }
                      }
                    j++;
                  }
                                                                    //  Now, is this actually what we want?
                if(strcmp(word0, "element") == 0 && strcmp(word1, "vertex") == 0)
                  (*info)[0] = (unsigned int)atoi(word2);
                if(strcmp(word0, "element") == 0 && strcmp(word1, "face") == 0)
                  (*info)[1] = (unsigned int)atoi(word2);
                if(strcmp(word0, "element") == 0 && strcmp(word1, "tetrahedra") == 0)
                  (*info)[2] = (unsigned int)atoi(word2);
              }

            linectr = 0;                                            //  Reset the line buffer's index, prepare for the next line

            if(strcmp(line, "end_header") == 0)                     //  We're done!
              {
                if(contentsLen > 0)
                  free(buffer);                                     //  Clean up
                (*info)[3] = ++i;                                   //  Point to the first byte of PLY data
                return true;                                        //  Let's get out of here!
              }
          }
        else
          {
            line[linectr] = buffer[i];                              //  Add this character to the line buffer
            linectr++;                                              //  Point to the next character
          }
        i++;
      }

    if(contentsLen > 0)
      free(buffer);                                                 //  Clean up

    return false;
  }

/* Compute and store the barycenter for this face. */
void computeFaceBarycenter(Vertex* V, Face** F, unsigned int f)
  {
    (*F)[f].barycenter[0] = (V[(*F)[f].v[0]].x + V[(*F)[f].v[1]].x + V[(*F)[f].v[2]].x) / 3.0;
    (*F)[f].barycenter[1] = (V[(*F)[f].v[0]].y + V[(*F)[f].v[1]].y + V[(*F)[f].v[2]].y) / 3.0;
    (*F)[f].barycenter[2] = (V[(*F)[f].v[0]].z + V[(*F)[f].v[1]].z + V[(*F)[f].v[2]].z) / 3.0;
    return;
  }

/* Area of a 3D triangle is 1/2 times the magnitude of the cross product of two of its vectors. */
void computeFaceArea(Vertex* V, Face** F, unsigned int f)
  {
    float* u;
    float* v;
    float* cross;

    if((u = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product input array u." << endl;
        exit(1);
      }
    if((v = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product input array v." << endl;
        exit(1);
      }
    if((cross = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product output array." << endl;
        exit(1);
      }

    u[0] = V[(*F)[f].v[0]].x - V[(*F)[f].v[1]].x;
    u[1] = V[(*F)[f].v[0]].y - V[(*F)[f].v[1]].y;
    u[2] = V[(*F)[f].v[0]].z - V[(*F)[f].v[1]].z;

    v[0] = V[(*F)[f].v[0]].x - V[(*F)[f].v[2]].x;
    v[1] = V[(*F)[f].v[0]].y - V[(*F)[f].v[2]].y;
    v[2] = V[(*F)[f].v[0]].z - V[(*F)[f].v[2]].z;

    crossProduct(u, v, &cross);

    (*F)[f].area = 0.5 * sqrt(pow(cross[0], 2) + pow(cross[1], 2) + pow(cross[2], 2));

    free(u);
    free(v);
    free(cross);

    return;
  }

/* Compute and store the unit normal for the Face at index 'f' in array 'F'. */
void computeFaceNormal(Vertex* V, Face** F, unsigned int f)
  {
    double Nx, Ny, Nz;
    double mag;

    //  (B - A) cross (C - A)

    //  [a]       [d]     [bf - ce]
    //  [b] cross [e]  =  [cd - af]
    //  [c]       [f]     [ae - bd]

    //             B                   A
    //  a = V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x
    //  b = V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y
    //  c = V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z

    //             C                   A
    //  d = V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x
    //  e = V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y
    //  f = V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z

    Nx = (V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y) * (V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z) - (V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z) * (V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y);
    Ny = (V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z) * (V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x) - (V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x) * (V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z);
    Nz = (V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x) * (V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y) - (V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y) * (V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x);
    mag = sqrt(pow(Nx, 2) + pow(Ny, 2) + pow(Nz, 2));

    (*F)[f].normal[0] = Nx / mag;
    (*F)[f].normal[1] = Ny / mag;
    (*F)[f].normal[2] = Nz / mag;

    return;
  }

/* Compute the cross-product, stick the result in the given array 'output' and return the magnitude. */
float crossProduct(float* u, float* v, float** output)
  {
    (*output)[0] = u[1] * v[2] - u[2] * v[1];
    (*output)[1] = u[2] * v[0] - u[0] * v[2];
    (*output)[2] = u[0] * v[1] - u[1] * v[0];
    return sqrt(pow((*output)[0], 2) + pow((*output)[1], 2) + pow((*output)[2], 2));
  }

/* Arrange both 'areas' and 'indices' according to values in 'areas'. */
void quicksort(double** areas, unsigned int** indices, unsigned int lo, unsigned int hi)
  {
    unsigned int p;

    if(lo < hi)
      {
        p = partition(false, areas, indices, lo, hi);               //  Sort ascending by largest area

        if(p > 0)                                                   //  PREVENT ROLL-OVER TO INT_MAX
          quicksort(areas, indices, lo, p - 1);                     //  Left side: start quicksort
        if(p < INT_MAX)                                             //  PREVENT ROLL-OVER TO 0
          quicksort(areas, indices, p + 1, hi);                     //  Right side: start quicksort
      }

    return;
  }

unsigned int partition(bool desc, double** areas, unsigned int** indices, unsigned int lo, unsigned int hi)
  {
    double pivot;
    unsigned int i = lo;
    unsigned int j;
    unsigned int tmpInt;
    double tmpDouble;
    bool trigger;

    pivot = (*areas)[hi];

    for(j = lo; j < hi; j++)
      {
        if(desc)
          trigger = ((*areas)[j] > pivot);                          //  SORT DESCENDING
        else
          trigger = ((*areas)[j] < pivot);                          //  SORT ASCENDING

        if(trigger)
          {
            tmpInt = (*indices)[i];                                 //  tmp gets [i]
            (*indices)[i] = (*indices)[j];                          //  [i] gets [j]
            (*indices)[j] = tmpInt;                                 //  [j] gets tmp

            tmpDouble = (*areas)[i];                                //  tmp gets [i]
            (*areas)[i] = (*areas)[j];                              //  [i] gets [j]
            (*areas)[j] = tmpDouble;                                //  [j] gets tmp

            i++;
          }
      }

    tmpInt = (*indices)[i];                                         //  tmp gets [i]
    (*indices)[i] = (*indices)[hi];                                 //  [i] gets [hi]
    (*indices)[hi] = tmpInt;                                        //  [hi] gets tmp

    tmpDouble = (*areas)[i];                                        //  tmp gets [i]
    (*areas)[i] = (*areas)[hi];                                     //  [i] gets [hi]
    (*areas)[hi] = tmpDouble;                                       //  [hi] gets tmp

    return i;
  }

void usage(void)
  {
    cout << "Usage:  ./areaTh ply-filename <area-threshold>" << endl;
    cout << " e.g.:  ./areaTh terrains.ply" << endl;
    cout << " e.g.:  ./areaTh terrains.ply 0.1" << endl;
    cout << "Writes a new file \"areath.ply\" equal to the original, minus the too-large triangles and any now-unnecessary vertices." << endl;
    return;
  }
