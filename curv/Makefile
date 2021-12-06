CC = gcc
_OBJS_RW = rw_mesh.o rw_mesh_header.o rw_mesh_mesh.o rw_mesh_off.o rw_mesh_string.o rw_mesh_vtk.o
OBJS_RW   = $(patsubst %,rw_mesh/%, $(_OBJS_RW))
_OBJS = curv.o dual.o main.o figures.o
OBJS = $(patsubst %,curvature/%, $(_OBJS))
_SRC_RW = rw_mesh.c rw_mesh_header.c rw_mesh_mesh.c rw_mesh_off.c rw_mesh_string.c rw_mesh_vtk.c 
_HDRS_RW  = rw_mesh.h rw_mesh_header.h rw_mesh_mesh.h rw_mesh_off.h rw_mesh_string.h rw_mesh_vtk.h
HDRS_RW   = $(patsubst %,rw_mesh/%, $(_HDRS_RW))
_HDRS = curvature.h
HDRS = $(patsubst %,curvature/%, $(_HDRS)) 
SRC_RW = $(patsubst %,rw_mesh/%, $(_SRC_RW))
_SRC = curv.c dual.c main.c figures.c
SRC = $(patsubst %,curvature/%, $(_SRC))
PROG = curv
LIB_RW = lib_rw_mesh.a
LIB = lib_curv.a

all: rw_mesh $(PROG)

main.o : $(HDRS) $(HDRS_RW) 

$(OBJS) : $(HDRS) 

$(OBJS_RW) : $(HDRS_RW)

rw_mesh: 
	cd rw_mesh && $(CC) -w -c $(_SRC_RW)
	ar rcv rw_mesh/librwmsh.a $(OBJS_RW)

$(PROG):   
	$(CC) $(SRC) -o $(PROG) -Lrw_mesh -lrwmsh -lm



