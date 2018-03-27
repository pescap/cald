import bempp.api 

_l_s = 0


def get_gmsh_file():
    """
    Return a 3-tuple (geo_file,geo_name,msh_name), where
    geo_file is a file descriptor to an empty .geo file, geo_name is
    the corresponding filename and msh_name is the name of the
    Gmsh .msh file that will be generated.

    """
    import os
    import tempfile

    geo, geo_name = tempfile.mkstemp(
        suffix='.geo', dir=bempp.api.TMP_PATH, text=True)
    geo_file = os.fdopen(geo, "w")
    msh_name = os.path.splitext(geo_name)[0] + ".msh"
    return (geo_file, geo_name, msh_name)


def __generate_grid_from_gmsh_string(gmsh_string):
    """Return a grid from a string containing a gmsh mesh"""

    import os
    import tempfile

    handle, fname = tempfile.mkstemp(
        suffix='.msh', dir=bempp.api.TMP_PATH, text=True)
    with os.fdopen(handle, "w") as f:
        f.write(gmsh_string)
    grid = bempp.api.import_grid(fname)
    os.remove(fname)
    return grid


def __generate_grid_from_geo_string(geo_string):
    """Helper routine that implements the grid generation
    """

    import os
    import subprocess

    def msh_from_string(geo_string):
        gmsh_command = bempp.api.GMSH_PATH
        if gmsh_command is None:
            raise RuntimeError("Gmsh is not found. Cannot generate mesh")
        f, geo_name, msh_name = get_gmsh_file()
        f.write(geo_string)
        f.close()

        fnull = open(os.devnull, 'w')
        cmd = gmsh_command + " -2 " + geo_name
        try:
            subprocess.check_call(
                cmd, shell=True, stdout=fnull, stderr=fnull)
        except:
            print("The following command failed: " + cmd)
            fnull.close()
            raise
        os.remove(geo_name)
        fnull.close()
        return msh_name

    msh_name = msh_from_string(geo_string)
    grid = bempp.api.import_grid(msh_name)
    os.remove(msh_name)
    return grid

def _make_word(geo,h):
    intro = "lc = "+str(h)+";\n"
    outro = "\nMesh.Algorithm = 6;"
    return __generate_grid_from_geo_string(intro+geo+outro)




def letters_m(h=0.1):
    return _make_word(_m(0,0),h)



_twelve = """Line(ls+1) = {ls+1,ls+2};
           Line(ls+2) = {ls+2,ls+3};
           Line(ls+3) = {ls+3,ls+4};
           Line(ls+4) = {ls+4,ls+5};
           Line(ls+5) = {ls+5,ls+6};
           Line(ls+6) = {ls+6,ls+7};
           Line(ls+7) = {ls+7,ls+8};
           Line(ls+8) = {ls+8,ls+9};
           Line(ls+9) = {ls+9,ls+10};
           Line(ls+10) = {ls+10,ls+11};
           Line(ls+11) = {ls+11,ls+12};
           Line(ls+12) = {ls+12,ls+1};

           Line(ls+13) = {ls+13,ls+14};
           Line(ls+14) = {ls+14,ls+15};
           Line(ls+15) = {ls+15,ls+16};
           Line(ls+16) = {ls+16,ls+17};
           Line(ls+17) = {ls+17,ls+18};
           Line(ls+18) = {ls+18,ls+19};
           Line(ls+19) = {ls+19,ls+20};
           Line(ls+20) = {ls+20,ls+21};
           Line(ls+21) = {ls+21,ls+22};
           Line(ls+22) = {ls+22,ls+23};
           Line(ls+23) = {ls+23,ls+24};
           Line(ls+24) = {ls+24,ls+13};

           Line(ls+25) = {ls+1,ls+13};
           Line(ls+26) = {ls+2,ls+14};
           Line(ls+27) = {ls+3,ls+15};
           Line(ls+28) = {ls+4,ls+16};
           Line(ls+29) = {ls+5,ls+17};
           Line(ls+30) = {ls+6,ls+18};
           Line(ls+31) = {ls+7,ls+19};
           Line(ls+32) = {ls+8,ls+20};
           Line(ls+33) = {ls+9,ls+21};
           Line(ls+34) = {ls+10,ls+22};
           Line(ls+35) = {ls+11,ls+23};
           Line(ls+36) = {ls+12,ls+24};


           Line Loop(ls+1) = {ls+1,ls+2,ls+3,ls+4,ls+5,ls+6,ls+7,ls+8,ls+9,ls+10,ls+11,ls+12};
           Line Loop(ls+2) = {-ls-24,-ls-23,-ls-22,-ls-21,-ls-20,-ls-19,-ls-18,-ls-17,-ls-16,-ls-15,-ls-14,-ls-13};

           Line Loop(ls+3) = {-ls-1,ls+25,ls+13,-ls-26};
           Line Loop(ls+4) = {-ls-2,ls+26,ls+14,-ls-27};
           Line Loop(ls+5) = {-ls-3,ls+27,ls+15,-ls-28};
           Line Loop(ls+6) = {-ls-4,ls+28,ls+16,-ls-29};
           Line Loop(ls+7) = {-ls-5,ls+29,ls+17,-ls-30};
           Line Loop(ls+8) = {-ls-6,ls+30,ls+18,-ls-31};
           Line Loop(ls+9) = {-ls-7,ls+31,ls+19,-ls-32};
           Line Loop(ls+10) = {-ls-8,ls+32,ls+20,-ls-33};
           Line Loop(ls+11) = {-ls-9,ls+33,ls+21,-ls-34};
           Line Loop(ls+12) = {-ls-10,ls+34,ls+22,-ls-35};
           Line Loop(ls+13) = {-ls-11,ls+35,ls+23,-ls-36};
           Line Loop(ls+14) = {-ls-12,ls+36,ls+24,-ls-25};

           Plane Surface(ls+1) = {ls+1};
           Plane Surface(ls+2) = {ls+2};
           Plane Surface(ls+3) = {ls+3};
           Plane Surface(ls+4) = {ls+4};
           Plane Surface(ls+5) = {ls+5};
           Plane Surface(ls+6) = {ls+6};
           Plane Surface(ls+7) = {ls+7};
           Plane Surface(ls+8) = {ls+8};
           Plane Surface(ls+9) = {ls+9};
           Plane Surface(ls+10) = {ls+10};
           Plane Surface(ls+11) = {ls+11};
           Plane Surface(ls+12) = {ls+12};
           Plane Surface(ls+14) = {ls+14};
           Plane Surface(ls+13) = {ls+13};"""


def _m(x=0,y=0):
    """Return gmsh code for a letter E with bottom left corner at (x,y)."""
    global _l_s
    out = "ls = "+str(_l_s)+""";
           Point(ls+1) = {"""+str(x)+","+str(y)+""",1,lc};
           Point(ls+2) = {"""+str(x+.5)+","+str(y)+""",1,lc};
           Point(ls+3) = {"""+str(x+.5)+","+str(y+3)+""",1,lc};
           Point(ls+4) = {"""+str(x+1.5)+","+str(y+1)+""",1,lc};
           Point(ls+5) = {"""+str(x+2.5)+","+str(y+3)+""",1,lc};
           Point(ls+6) = {"""+str(x+2.5)+","+str(y)+""",1,lc};
           Point(ls+7) = {"""+str(x+3)+","+str(y)+""",1,lc};
           Point(ls+8) = {"""+str(x+3)+","+str(y+3.92)+""",1,lc};
           Point(ls+9) = {"""+str(x+2.3)+","+str(y+3.92)+""",1,lc};
           Point(ls+10) = {"""+str(x+1.5)+","+str(y+2.32)+""",1,lc};
           Point(ls+11) = {"""+str(x+0.9)+","+str(y+3.92)+""",1,lc};
           Point(ls+12) = {"""+str(x)+","+str(y+3.92)+""",1,lc};

           Point(ls+13) = {"""+str(x)+","+str(y)+""",0,lc};
           Point(ls+14) = {"""+str(x+.5)+","+str(y)+""",0,lc};
           Point(ls+15) = {"""+str(x+.5)+","+str(y+3)+""",0,lc};
           Point(ls+16) = {"""+str(x+1.5)+","+str(y+1)+""",0,lc};
           Point(ls+17) = {"""+str(x+2.5)+","+str(y+3)+""",0,lc};
           Point(ls+18) = {"""+str(x+2.5)+","+str(y)+""",0,lc};
           Point(ls+19) = {"""+str(x+3)+","+str(y)+""",0,lc};
           Point(ls+20) = {"""+str(x+3)+","+str(y+3.92)+""",0,lc};
           Point(ls+21) = {"""+str(x+2.3)+","+str(y+3.92)+""",0,lc};
           Point(ls+22) = {"""+str(x+1.5)+","+str(y+2.32)+""",0,lc};
           Point(ls+23) = {"""+str(x+0.9)+","+str(y+3.92)+""",0,lc};
           Point(ls+24) = {"""+str(x)+","+str(y+3.92)+""",0,lc};
           """ + _twelve

    _l_s += 36
    return out


def cube_mentor(lx=1, ly=1,lz=1, origin=(0, 0, 0), h=0.1):
    cube_stub = """
    Point(1) = {orig0,orig1,orig2,cl};
    Point(2) = {orig0+lx,orig1,orig2,cl};
    Point(3) = {orig0+lx,orig1+ly,orig2,cl};
    Point(4) = {orig0,orig1+ly,orig2,cl};
    Point(5) = {orig0,orig1,orig2+lz,cl};
    Point(6) = {orig0+lx,orig1,orig2+lz,cl};
    Point(7) = {orig0+lx,orig1+ly,orig2+lz,cl};
    Point(8) = {orig0,orig1+ly,orig2+lz,cl};

    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Line(5) = {1,5};
    Line(6) = {2,6};
    Line(7) = {3,7};
    Line(8) = {4,8};
    Line(9) = {5,6};
    Line(10) = {6,7};
    Line(11) = {7,8};
    Line(12) = {8,5};

    Line Loop(1) = {-1,-4,-3,-2};
    Line Loop(2) = {1,6,-9,-5};
    Line Loop(3) = {2,7,-10,-6};
    Line Loop(4) = {3,8,-11,-7};
    Line Loop(5) = {4,5,-12,-8};
    Line Loop(6) = {9,10,11,12};

    Plane Surface(1) = {1};
    Plane Surface(2) = {2};
    Plane Surface(3) = {3};
    Plane Surface(4) = {4};
    Plane Surface(5) = {5};
    Plane Surface(6) = {6};

    Physical Surface(1) = {1};
    Physical Surface(2) = {2};
    Physical Surface(3) = {3};
    Physical Surface(4) = {4};
    Physical Surface(5) = {5};
    Physical Surface(6) = {6};

    Surface Loop (1) = {1,2,3,4,5,6};

    Volume (1) = {1};

    Mesh.Algorithm = 6;
    """

    cube_geometry = (
        "lx = " + str(lx) + ";\n" +
        "ly = " + str(ly) + ";\n" +
        "lz = " + str(lz) + ";\n" +
        "orig0 = " + str(origin[0]) + ";\n" +
        "orig1 = " + str(origin[1]) + ";\n" +
        "orig2 = " + str(origin[2]) + ";\n" +
        "cl = " + str(h) + ";\n" + cube_stub)

    return __generate_grid_from_geo_string(cube_geometry)




def cube_mentor_graded(lx=2, ly=2,lz=0.1, origin=(-1, -1, 0), h=0.1):
    cube_stub = """
    Point(1) = {orig0,orig1,orig2,cl};
    Point(2) = {orig0+lx,orig1,orig2,cl};
    Point(3) = {orig0+lx,orig1+ly,orig2,cl};
    Point(4) = {orig0,orig1+ly,orig2,cl};
    Point(5) = {orig0,orig1,orig2+lz,cl};
    Point(6) = {orig0+lx,orig1,orig2+lz,cl};
    Point(7) = {orig0+lx,orig1+ly,orig2+lz,cl};
    Point(8) = {orig0,orig1+ly,orig2+lz,cl};

    Line(1) = {1,2};
    Line(2) = {2,3};
    Line(3) = {3,4};
    Line(4) = {4,1};
    Line(5) = {1,5};
    Line(6) = {2,6};
    Line(7) = {3,7};
    Line(8) = {4,8};
    Line(9) = {5,6};
    Line(10) = {6,7};
    Line(11) = {7,8};
    Line(12) = {8,5};

    Line Loop(1) = {-1,-4,-3,-2};
    Line Loop(2) = {1,6,-9,-5};
    Line Loop(3) = {2,7,-10,-6};
    Line Loop(4) = {3,8,-11,-7};
    Line Loop(5) = {4,5,-12,-8};
    Line Loop(6) = {9,10,11,12};

    Plane Surface(1) = {1};
    Plane Surface(2) = {2};
    Plane Surface(3) = {3};
    Plane Surface(4) = {4};
    Plane Surface(5) = {5};
    Plane Surface(6) = {6};

    Physical Surface(1) = {1};
    Physical Surface(2) = {2};
    Physical Surface(3) = {3};
    Physical Surface(4) = {4};
    Physical Surface(5) = {5};
    Physical Surface(6) = {6};

    Surface Loop (1) = {1,2,3,4,5,6};

    Field[1] = Attractor;
    Field[1].NNodesByEdge = 100;
    Field[1].EdgesList = {1,5,9};

    Field[2] = Attractor;
    Field[2].NodesList = {1,2,3,4,5,6,7,8};
    

    Field[3] = Threshold;
    Field[3].IField = 1;
    Field[3].LcMin = cl / 5;
    Field[3].LcMax = 2*cl;
    Field[3].DistMin = 0.01;
    Field[3].DistMax = 0.05;

    Field[4] = Threshold;
    Field[4].IField = 2;
    Field[4].LcMin = cl / 10;
    Field[4].LcMax = 2*cl;
    Field[4].DistMin = 0.1;
    Field[4].DistMax = 0.5;

    Field[5] = Min;
    Field[5].FieldsList = {3,4};

    Background Field = 5;
    Mesh.Algorithm = 3;


    """

    cube_geometry = (
        "lx = " + str(lx) + ";\n" +
        "ly = " + str(ly) + ";\n" +
        "lz = " + str(lz) + ";\n" +
        "orig0 = " + str(origin[0]) + ";\n" +
        "orig1 = " + str(origin[1]) + ";\n" +
        "orig2 = " + str(origin[2]) + ";\n" +
        "cl = " + str(h) + ";\n" + cube_stub)

    return __generate_grid_from_geo_string(cube_geometry)


def destroyer(origin=(0, 0, 0), h=100):
    cube_stub = """
      Point(1) = {39.8035998778,-24.7302271217,-41.0004006351,lc};
      Point(2) = {-66.3077627774,-23.8877610579,-41.0004006351,lc};
      Point(3) = {-136.679671026,-22.057582714,-41.0004006351,lc};
      Point(4) = {225.839375183,-21.306256563,-41.0004006351,lc};
      Point(5) = {313.903782423,-18.1035498912,-41.0004006351,lc};
      Point(6) = {-239.436261148,-13.3910868899,-41.0004006351,lc};
      Point(7) = {350.075915239,-12.9798537372,-41.0004006351,lc};
      Point(8) = {-332.894256544,-2.2486622129e-06,-41.0004020223,lc};
      Point(9) = {359.148439312,-2.2486622129e-06,-41.0004020223,lc};
      Point(10) = {-239.436261148,13.3910810054,-41.0004034096,lc};
      Point(11) = {350.075915239,12.9798492398,-41.0004034096,lc};
      Point(12) = {313.903782423,18.1035453938,-41.0004034096,lc};
      Point(13) = {225.839375183,21.3062520657,-41.0004034096,lc};
      Point(14) = {-136.679671026,22.057579604,-41.0004034096,lc};
      Point(15) = {-66.3077627774,23.8877565606,-41.0004034096,lc};
      Point(16) = {39.8035998778,24.7302240116,-41.0004034096,lc};
      Point(17) = {241.123930737,-23.8052206111,-27.5030274813,lc};
      Point(18) = {309.882187123,-20.6518980966,-27.5030274813,lc};
      Point(19) = {344.822475213,-14.8069246759,-27.5030274813,lc};
      Point(20) = {352.886453521,-8.61427812973e-07,-27.5030288686,lc};
      Point(21) = {344.822475213,14.8069229531,-27.5030302558,lc};
      Point(22) = {309.882187123,20.6518949865,-27.5030302558,lc};
      Point(23) = {241.123930737,23.805217501,-27.5030302558,lc};
      Point(24) = {39.8035998778,-28.2113531249,-18.058501244,lc};
      Point(25) = {-66.3086819866,-27.3661932652,-18.0762062516,lc};
      Point(26) = {225.839375183,-24.3053186123,-18.058501244,lc};
      Point(27) = {84.2010671707,-21.6785888888,-18.0586385802,lc};
      Point(28) = {185.240564911,-21.6785888888,-18.0586385802,lc};
      Point(29) = {200.143639329,-7.69469875861,-18.058501244,lc};
      Point(30) = {213.505465809,-0.0194068834484,-18.058501244,lc};
      Point(31) = {200.178311865,7.71410444501,-18.058501244,lc};
      Point(32) = {84.2010671707,21.6785885531,-18.0586413547,lc};
      Point(33) = {185.240564911,21.6785885531,-18.0586413547,lc};
      Point(34) = {225.839375183,24.3053155022,-18.0585040185,lc};
      Point(35) = {-66.3086819773,27.3661887992,-18.0762088203,lc};
      Point(36) = {39.8035998778,28.2113500148,-18.0585040185,lc};
      Point(37) = {-136.679671026,-25.1624893607,-17.0806397155,lc};
      Point(38) = {-136.679671026,25.1624862506,-17.0806411027,lc};
      Point(39) = {-239.436261148,-15.27608735,-14.8588437133,lc};
      Point(40) = {-239.436261148,15.2760842399,-14.8588464878,lc};
      Point(41) = {-228.704626928,-6.85997068825,-14.3558339072,lc};
      Point(42) = {-214.984878713,-6.85997068825,-14.3558339072,lc};
      Point(43) = {-228.704626928,6.85997035263,-14.3558339072,lc};
      Point(44) = {-214.984878713,6.85997035263,-14.3558339072,lc};
      Point(45) = {-340.152259988,-8.61427812973e-07,-12.2544498509,lc};
      Point(46) = {204.609138539,-2.2486622129e-06,-10.3993026751,lc};
      Point(47) = {-228.704626928,-0.368109058386,-7.56601513671,lc};
      Point(48) = {-228.704626928,0.368107335531,-7.56601513671,lc};
      Point(49) = {-240.15902553,-0.265761678828,-7.3887265804,lc};
      Point(50) = {-240.15902553,0.265759955973,-7.3887265804,lc};
      Point(51) = {-228.704626928,-0.736218642579,-6.92844081927,lc};
      Point(52) = {-240.15902553,-0.531523883464,-6.92844081927,lc};
      Point(53) = {-240.15902553,0.531522160608,-6.92844081927,lc};
      Point(54) = {-228.704626928,0.736216919724,-6.92844081927,lc};
      Point(55) = {-240.15902553,-0.265761678828,-6.46815644537,lc};
      Point(56) = {-240.15902553,0.265759955973,-6.46815644537,lc};
      Point(57) = {-228.704626928,-0.368109058386,-6.29086788906,lc};
      Point(58) = {-228.704626928,0.368107335531,-6.29087066353,lc};
      Point(59) = {30.9168375885,-23.722814726,-4.58124160183,lc};
      Point(60) = {84.2010671707,-21.6785888888,-4.58387734719,lc};
      Point(61) = {185.240564911,-21.6785888888,-4.58387734719,lc};
      Point(62) = {-2.65378958757,-23.3187133453,-4.54850286999,lc};
      Point(63) = {-88.2640779061,-22.9147506881,-4.51312700555,lc};
      Point(64) = {-131.08217606,-21.159069606,-4.51312700555,lc};
      Point(65) = {94.250886781,-17.4264380061,-4.58387734719,lc};
      Point(66) = {175.190047522,-17.4264380061,-4.58387734719,lc};
      Point(67) = {11.7114928961,-7.69469737137,-4.58124160183,lc};
      Point(68) = {25.0726965084,-0.019405496214,-4.58124437629,lc};
      Point(69) = {11.7451055857,7.71410444501,-4.58124437629,lc};
      Point(70) = {94.250886781,17.4264376705,-4.58387873442,lc};
      Point(71) = {175.190047522,17.4264376705,-4.58387873442,lc};
      Point(72) = {84.2010671707,21.6785885531,-4.58387873442,lc};
      Point(73) = {185.240564911,21.6785885531,-4.58387873442,lc};
      Point(74) = {30.9168375885,23.7228171649,-4.58124437629,lc};
      Point(75) = {-131.08217606,21.1590678831,-4.51312956745,lc};
      Point(76) = {-2.65378958757,23.3187157842,-4.54850564446,lc};
      Point(77) = {-228.704626928,-5.45736573117,-4.0113643231,lc};
      Point(78) = {-215.62993716,-5.45736573117,-4.0113643231,lc};
      Point(79) = {-88.2640779061,22.9147517397,-4.51312978002,lc};
      Point(80) = {-228.704626928,5.45736539555,-4.01136709757,lc};
      Point(81) = {-215.62993716,5.45736539555,-4.01136709757,lc};
      Point(82) = {16.1764455358,-8.61427812973e-07,3.05423387159,lc};
      Point(83) = {98.6319117394,-13.0249208211,11.324925364,lc};
      Point(84) = {159.128643478,-13.0249208211,11.324925364,lc};
      Point(85) = {100.6686692893,-10.9906261949,11.324925364,lc};
      Point(86) = {149.716263623,-10.9906261949,11.324925364,lc};
      Point(87) = {100.6685305658,10.9906286338,11.324919815,lc};
      Point(88) = {149.716263623,10.9906286338,11.324919815,lc};
      Point(89) = {98.6319117394,13.02492326,11.324919815,lc};
      Point(90) = {159.128643478,13.02492326,11.324919815,lc};
      Point(91) = {-13.130906525,-14.8606078728,23.7799347513,lc};
      Point(92) = {-77.7868499899,-14.6031371681,23.8046261364,lc};
      Point(93) = {102.081401862,-7.93484013111,23.7803523089,lc};
      Point(94) = {138.936059552,-7.93484013111,23.7803523089,lc};
      Point(95) = {-47.4226500485,-6.09743429606,23.7799292024,lc};
      Point(96) = {-19.102118276,-6.09743429606,23.7799292024,lc};
      Point(97) = {-66.6927217106,-4.28175255171,23.8046219747,lc};
      Point(98) = {-52.9728361591,-1.19593030124e-05,23.8046219747,lc};
      Point(99) = {-66.6927217106,4.2817563778,23.8046219747,lc};
      Point(100) = {-47.4226500485,6.09743673491,23.7799292024,lc};
      Point(101) = {-19.102118276,6.09743673491,23.7799292024,lc};
      Point(102) = {102.081401862,7.93484256996,23.7803467599,lc};
      Point(103) = {138.936059552,7.93484256996,23.7803467599,lc};
      Point(104) = {-13.130906525,14.8606116988,23.7799292024,lc};
      Point(105) = {-77.7868499899,14.603139607,23.8046219747,lc};
      Point(106) = {-45.5292096543,-4.4021228806,30.4966421072,lc};
      Point(107) = {-25.0829019444,-4.4021228806,30.4966421072,lc};
      Point(108) = {-45.5292096543,4.40212670668,30.4966421072,lc};
      Point(109) = {-25.0829019444,4.40212670668,30.4966421072,lc};
      Point(110) = {-56.3841828848,-23.8878956196,42.2137787427,lc};
      Point(111) = {-56.3841828848,23.8879022202,42.2137787427,lc};
      Point(112) = {-57.7603194095,-23.8878956196,44.5974622249,lc};
      Point(113) = {-55.0079076366,-23.8878956196,44.5974622249,lc};
      Point(114) = {-57.7603194095,23.8879022202,44.5974622249,lc};
      Point(115) = {-55.0079076366,23.8879022202,44.5974622249,lc};
      Point(116) = {-56.3841828848,-23.8878956196,58.5595616771,lc};
      Point(117) = {-56.3841828848,23.8879022202,58.5595519664,lc};
      Point(118) = {-57.7603194095,-23.8878956196,60.9432465465,lc};
      Point(119) = {-55.0079076366,-23.8878956196,60.9432465465,lc};
      Point(120) = {-57.7603194095,23.8879022202,60.9432368358,lc};
      Point(121) = {-55.0079076366,23.8879022202,60.9432368358,lc};
      Point(122) = {-66.6927217106,-4.2817563778,72.3282695561,lc};
      Point(123) = {-52.9728361591,-1.19593030124e-05,72.3282695561,lc};
      Point(124) = {-66.6927217106,4.2817563778,72.3282695561,lc};
      Point(125) = {-55.0079076366,0.63510300454945878,60.9432465465,lc};
      Point(126) = {-55.0079076366,-0.63510300454945878,60.9432465465,lc};
      Point(127) = {-57.7603194095,1.4940890023287405,60.9432465465,lc};
      Point(128) = {-56.3841828848,1.064617650190002,58.5595616771,lc};
      Point(129) = {-56.3841828848,-1.064617650190002,58.5595616771,lc};
      Point(130) = {-57.7603194095,-1.4940890023287405,60.9432465465,lc};
      Point(131) = {-57.7603194095,-1.4940890023287405,44.5974622249,lc};
      Point(132) = {-55.0079076366,-0.63510300454945878,44.5974622249,lc};
      Point(133) = {-56.3841828848,-1.064617650190002,42.2137787427,lc};
      Point(134) = {-57.7603194095,1.4940890023287405,44.5974622249,lc};
      Point(135) = {-55.0079076366,0.63510300454945878,44.5974622249,lc};
      Point(136) = {-56.3841828848,1.064617650190002,42.2137787427,lc};
      Line(1) = {1,2};
      Line(2) = {4,1};
      Line(3) = {2,3};
      Line(5) = {16,1};
      Line(6) = {3,6};
      Line(7) = {5,4};
      Line(8) = {24,1};
      Line(9) = {15,2};
      Line(11) = {7,5};
      Line(12) = {14,3};
      Line(14) = {6,8};
      Line(17) = {2,25};
      Line(18) = {13,4};
      Line(20) = {10,6};
      Line(21) = {12,5};
      Line(22) = {9,7};
      Line(23) = {4,17};
      Line(26) = {8,10};
      Line(29) = {18,5};
      Line(31) = {11,9};
      Line(32) = {26,4};
      Line(33) = {3,37};
      Line(35) = {12,11};
      Line(36) = {7,19};
      Line(37) = {10,14};
      Line(39) = {13,12};
      Line(40) = {20,9};
      Line(41) = {16,13};
      Line(42) = {14,15};
      Line(45) = {11,21};
      Line(46) = {6,39};
      Line(47) = {15,16};
      Line(49) = {12,22};
      Line(51) = {13,23};
      Line(52) = {17,18};
      Line(55) = {18,19};
      Line(56) = {8,45};
      Line(59) = {19,20};
      Line(61) = {17,23};
      Line(62) = {18,22};
      Line(64) = {10,40};
      Line(65) = {20,21};
      Line(66) = {26,17};
      Line(67) = {34,13};
      Line(68) = {21,22};
      Line(71) = {22,23};
      Line(72) = {15,35};
      Line(73) = {14,38};
      Line(75) = {36,16};
      Line(77) = {24,25};
      Line(78) = {24,26};
      Line(79) = {28,27};
      Line(80) = {23,34};
      Line(82) = {27,32};
      Line(83) = {36,24};
      Line(84) = {29,30};
      Line(85) = {26,34};
      Line(87) = {31,29};
      Line(88) = {33,28};
      Line(89) = {25,37};
      Line(90) = {30,31};
      Line(91) = {33,32};
      Line(92) = {34,36};
      Line(93) = {35,36};
      Line(94) = {38,35};
      Line(95) = {46,29};
      Line(96) = {30,46};
      Line(97) = {37,38};
      Line(98) = {24,59};
      Line(99) = {31,46};
      Line(100) = {39,37};
      Line(102) = {24,62};
      Line(103) = {38,40};
      Line(104) = {62,25};
      Line(105) = {39,40};
      Line(106) = {25,63};
      Line(107) = {60,27};
      Line(109) = {28,61};
      Line(110) = {42,41};
      Line(111) = {45,39};
      Line(112) = {43,41};
      Line(113) = {40,45};
      Line(115) = {44,42};
      Line(116) = {43,44};
      Line(120) = {47,48};
      Line(121) = {49,47};
      Line(122) = {32,72};
      Line(123) = {63,37};
      Line(126) = {64,37};
      Line(127) = {51,47};
      Line(128) = {50,48};
      Line(129) = {73,33};
      Line(131) = {49,50};
      Line(133) = {52,49};
      Line(134) = {48,54};
      Line(135) = {50,53};
      Line(136) = {52,51};
      Line(137) = {35,76};
      Line(138) = {74,36};
      Line(142) = {76,36};
      Line(144) = {35,79};
      Line(148) = {38,75};
      Line(149) = {55,52};
      Line(150) = {54,53};
      Line(151) = {57,51};
      //Line(152) = {52,57};
      Line(153) = {53,56};
      Line(154) = {79,38};
      Line(156) = {56,55};
      Line(157) = {54,58};
      Line(158) = {55,57};
      Line(159) = {41,77};
      Line(162) = {56,58};
      Line(163) = {42,78};
      Line(165) = {58,57};
      Line(167) = {80,43};
      Line(169) = {81,44};
      Line(170) = {62,59};
      Line(171) = {60,61};
      Line(172) = {63,62};
      Line(173) = {64,63};
      Line(174) = {66,65};
      Line(175) = {72,60};
      Line(176) = {59,74};
      Line(178) = {61,73};
      Line(179) = {65,70};
      Line(180) = {67,68};
      Line(183) = {69,67};
      Line(184) = {71,66};
      Line(185) = {68,69};
      Line(186) = {76,62};
      Line(187) = {75,64};
      Line(188) = {71,70};
      Line(189) = {79,63};
      Line(191) = {72,73};
      Line(192) = {83,65};
      Line(194) = {82,67};
      Line(195) = {66,84};
      Line(196) = {68,82};
      Line(197) = {74,76};
      Line(198) = {62,91};
      Line(199) = {69,82};
      Line(201) = {63,92};
      Line(203) = {75,79};
      Line(205) = {79,76};
      Line(206) = {78,77};
      Line(207) = {77,80};
      Line(208) = {70,89};
      Line(211) = {78,81};
      Line(212) = {90,71};
      Line(213) = {80,81};
      Line(216) = {83,84};
      Line(217) = {86,85};
      Line(218) = {89,83};
      Line(219) = {85,87};
      Line(222) = {84,90};
      Line(223) = {88,86};
      Line(224) = {88,87};
      Line(225) = {104,76};
      Line(226) = {93,85};
      Line(228) = {89,90};
      Line(229) = {86,94};
      Line(233) = {105,79};
      Line(234) = {92,91};
      Line(235) = {93,94};
      Line(236) = {87,102};
      Line(238) = {103,88};
      Line(239) = {96,95};
      Line(240) = {91,104};
      Line(241) = {102,93};
      Line(242) = {95,100};
      Line(243) = {98,97};
      Line(247) = {97,99};
      Line(248) = {92,105};
      Line(249) = {94,103};
      Line(250) = {96,101};
      Line(251) = {98,99};
      Line(252) = {106,95};
      Line(253) = {101,100};
      Line(255) = {96,107};
      Line(256) = {102,103};
      Line(258) = {100,108};
      Line(261) = {105,104};
      Line(262) = {101,109};
      Line(263) = {106,107};
      Line(264) = {108,106};
      Line(266) = {109,107};
      Line(267) = {108,109};
      Line(268) = {97,122};
      Line(271) = {98,123};
      Line(273) = {111,136};
      Line(274) = {99,124};
      Line(275) = {110,112};
      Line(276) = {113,110};
      Line(279) = {111,114};
      Line(280) = {112,113};
      Line(281) = {115,111};
      Line(282) = {114,134};
      Line(284) = {115,135};
      Line(285) = {114,115};
      Line(286) = {117,128};
      Line(287) = {116,118};
      Line(288) = {119,116};
      Line(291) = {117,120};
      Line(292) = {118,119};
      Line(293) = {121,117};
      Line(294) = {120,127};
      Line(296) = {121,125};
      Line(297) = {120,121};
      Line(298) = {123,122};
      Line(299) = {124,122};
      Line(300) = {123,124};
      Line(301) = {130,126};
      Line(302) = {126,129};
      Line(303) = {129,130};
      Line(304) = {127,125};
      Line(305) = {125,128};
      Line(306) = {128,127};
      Line(307) = {126,119};
      Line(308) = {130,118};
      Line(309) = {129,116};
      Line(310) = {131,132};
      Line(311) = {132,133};
      Line(312) = {133,131};
      Line(313) = {134,135};
      Line(314) = {135,136};
      Line(315) = {136,134};
      Line(316) = {131,112};
      Line(317) = {132,113};
      Line(318) = {133,110};
      Line Loop(500000) = {201,248,233,189};
      Line Loop(500006) = {-253,-250,239,242};
      Line Loop(500008) = {263,-266,-267,264};
      Line Loop(500010) = {252,-239,255,-263};
      Line Loop(500016) = {262,266,-255,250};
      Line Loop(500022) = {-264,-258,-242,-252};
      Line Loop(500024) = {-262,253,258,267};
      Line Loop(500032) = {163,206,-159,-110};
      Line Loop(500034) = {-116,-167,213,169};
      Line Loop(500036) = {211,-213,-207,-206};
      Line Loop(500038) = {115,110,-112,116};
      Line Loop(500040) = {-169,-211,-163,-115};
      Line Loop(500042) = {159,207,167,112};
      Line Loop(500058) = {85,92,83,78};
      Line Loop(500060) = {66,61,80,-85};
      Line Loop(500062) = {52,62,71,-61};
      Line Loop(500064) = {103,-105,100,97};
      Line Loop(500066) = {113,111,105};
      Line Loop(500076) = {59,65,68,-62,55};
      Line Loop(500084) = {18,2,-5,41};
      Line Loop(500086) = {1,-9,47,5};
      Line Loop(500088) = {-18,39,21,7};
      Line Loop(500090) = {3,-12,42,9};
      Line Loop(500092) = {6,-20,37,12};
      Line Loop(500094) = {31,22,11,-21,35};
      Line Loop(500096) = {14,26,20};
      Line Loop(500114) = {67,-41,-75,-92};
      Line Loop(500116) = {45,-65,40,-31};
      Line Loop(500118) = {72,93,75,-47};
      Line Loop(500120) = {-39,51,-71,-49};
      Line Loop(500122) = {73,94,-72,-42};
      Line Loop(500124) = {64,-103,-73,-37};
      Line Loop(500126) = {-35,49,-68,-45};
      Line Loop(500128) = {-113,-64,-26,56};
      Line Loop(500144) = {-67,-80,-51};
      Line Loop(500150) = {-2,-32,-78,8};
      Line Loop(500152) = {-14,46,-111,-56};
      Line Loop(500154) = {77,-17,-1,-8};
      Line Loop(500156) = {23,-66,32};
      Line Loop(500158) = {-3,17,89,-33};
      Line Loop(500160) = {-6,33,-100,-46};
      Line Loop(500162) = {-11,36,-55,29};
      Line Loop(500164) = {-22,-40,-59,-36};
      Line Loop(500180) = {-7,-29,-52,-23};
      Line Loop(500186) = {133,121,-127,-136};
      Line Loop(500188) = {131,128,-120,-121};
      Line Loop(500190) = {135,-150,-134,-128};
      Line Loop(500192) = {153,162,-157,150};
      Line Loop(500194) = {156,158,-165,-162};
      //Line Loop(500196) = {149,152,-158};
      Line Loop(500196) = {136,-151,-158,149};
      Line Loop(500198) = {134,157,165,151,127,120};
      Line Loop(500200) = {-153,-135,-131,-133,-149,-156};
      //Line Loop(500214) = {136,-151,-152};
      Line Loop(500228) = {271,298,-268,-243};
      Line Loop(500230) = {-300,-271,251,274};
      Line Loop(500232) = {299,-298,300};
      Line Loop(500234) = {247,-251,243};
      Line Loop(500236) = {268,-299,-274,-247};
      Line Loop(500244) = {284,314,-273,-281};
      Line Loop(500246) = {282,313,-284,-285};
      Line Loop(500248) = {273,315,-282,-279};
      Line Loop(500250) = {-276,-280,-275};
      Line Loop(500252) = {279,285,281};
      Line Loop(500260) = {296,305,-286,-293};
      Line Loop(500262) = {294,304,-296,-297};
      Line Loop(500264) = {286,306,-294,-291};
      Line Loop(500266) = {-288,-292,-287};
      Line Loop(500268) = {291,297,293};
      //Line Loop(500270) = {-294,-291,-289};
      Line Loop(500278) = {197,186,170,176};
      Line Loop(500280) = {-203,187,173,-189};
      Line Loop(500282) = {137,142,-93};
      Line Loop(500284) = {102,104,-77};
      Line Loop(500286) = {106,123,-89};
      Line Loop(500288) = {148,203,154};
      Line Loop(500290) = {-261,-248,234,240};
      Line Loop(500294) = {198,-234,-201,172};
      Line Loop(500296) = {-172,-106,-104};
      Line Loop(500298) = {144,205,-137};
      Line Loop(500308) = {225,-205,-233,261};
      Line Loop(500310) = {-144,-94,-154};
      Line Loop(500312) = {-197,138,-142};
      Line Loop(500314) = {98,-170,-102};
      Line Loop(500316) = {-173,126,-123};
      Line Loop(500318) = {-98,-83,-138,-176};
      Line Loop(500320) = {-126,-187,-148,-97};
      Line Loop(500322) = {-198,-186,-225,-240};
      Line Loop(500332) = {84,96,95};
      Line Loop(500334) = {90,99,-96};
      Line Loop(500336) = {87,-95,-99};
      Line Loop(500338) = {-87,-90,-84};
      Line Loop(500340) = {-91,88,79,82};
      Line Loop(500342) = {171,178,-191,175};
      Line Loop(500344) = {107,-79,109,-171};
      Line Loop(500352) = {129,91,122,191};
      Line Loop(500356) = {-188,184,174,179};
      Line Loop(500358) = {216,222,-228,218};
      Line Loop(500360) = {192,-174,195,-216};
      Line Loop(500368) = {212,188,208,228};
      Line Loop(500372) = {-224,223,217,219};
      Line Loop(500374) = {235,249,-256,241};
      Line Loop(500376) = {226,-217,229,-235};
      Line Loop(500378) = {236,256,238,224};
      Line Loop(500388) = {180,196,194};
      Line Loop(500390) = {185,199,-196};
      Line Loop(500392) = {183,-194,-199};
      Line Loop(500394) = {-183,-185,-180};
      Line Loop(500396) = {-107,117,-82};
      Line Loop(500398) = {-129,-178,-109,-88};
      Line Loop(500400) = {-175,-122,-82,-107};
      Line Loop(500406) = {-212,-222,-195,-184};
      Line Loop(500408) = {-218,-208,-179,-192};
      Line Loop(500414) = {-238,-249,-229,-223};
      Line Loop(500416) = {-241,-236,-219,-226};
      Line Loop(500420) = {303,301,302};
      Line Loop(500422) = {-308,-303,309,287};
      Line Loop(500424) = {-309,-302,307,288};
      Line Loop(500426) = {308,292,-307,-301};
      Line Loop(500428) = {305,306,304};
      Line Loop(500430) = {273,315,-282,-279};
      Line Loop(500432) = {282,313,-284,-285};
      Line Loop(500434) = {284,314,-273,-281};
      Line Loop(500436) = {-316,-312,318,275};
      Line Loop(500438) = {-318,-311,317,276};
      Line Loop(500440) = {316,280,-317,-310};
      Line Loop(500442) = {310,311,312};
      Line Loop(500444) = {313,314,315};
      Plane Surface(500001) = {500000};
      Plane Surface(500009) = {500008};
      Plane Surface(500011) = {500010};
      Plane Surface(500017) = {500016};
      Plane Surface(500023) = {500022};
      Plane Surface(500025) = {500024};
      Plane Surface(500033) = {500032};
      Plane Surface(500035) = {500034};
      Plane Surface(500037) = {500036};
      Plane Surface(500041) = {500040};
      Plane Surface(500043) = {500042,500198};
      Plane Surface(500059) = {500058,500338, 500340};
      Plane Surface(500061) = {500060};
      Plane Surface(500063) = {500062};
      Plane Surface(500065) = {500064,500038};
      Plane Surface(500067) = {500066};
      Plane Surface(500077) = {500076};
      Plane Surface(500085) = {500084};
      Plane Surface(500087) = {500086};
      Plane Surface(500089) = {500088};
      Plane Surface(500091) = {500090};
      Plane Surface(500093) = {500092};
      Plane Surface(500095) = {500094};
      Plane Surface(500097) = {500096};
      Plane Surface(500115) = {500114};
      Plane Surface(500117) = {500116};
      Plane Surface(500119) = {500118};
      Plane Surface(500121) = {500120};
      Plane Surface(500123) = {500122};
      Plane Surface(500125) = {500124};
      Plane Surface(500127) = {500126};
      Plane Surface(500129) = {500128};
      Plane Surface(500145) = {500144};
      Plane Surface(500151) = {500150};
      Plane Surface(500153) = {500152};
      Plane Surface(500155) = {500154};
      Plane Surface(500157) = {500156};
      Plane Surface(500159) = {500158};
      Plane Surface(500161) = {500160};
      Plane Surface(500163) = {500162};
      Plane Surface(500165) = {500164};
      Plane Surface(500181) = {500180};
      Plane Surface(500187) = {500186};
      Plane Surface(500189) = {500188};
      Plane Surface(500191) = {500190};
      Plane Surface(500193) = {500192};
      Plane Surface(500195) = {500194};
      Plane Surface(500197) = {500196};
      //Plane Surface(500199) = {500198};
      //Plane Surface(500201) = {500200};
      Plane Surface(500229) = {500228,500420,500442};
      Plane Surface(500231) = {500230,500444,500428};
      Plane Surface(500233) = {500232};
      Plane Surface(500237) = {500236};
      Plane Surface(500245) = {500244};
      Plane Surface(500247) = {500246};
      Plane Surface(500249) = {500248};
      Plane Surface(500251) = {500250};
      Plane Surface(500253) = {500252};
      //Plane Surface(500255) = {500254};
      Plane Surface(500261) = {500260};
      Plane Surface(500263) = {500262};
      Plane Surface(500265) = {500264};
      Plane Surface(500267) = {500266};
      Plane Surface(500269) = {500268};
      //Plane Surface(500271) = {500270};
      Plane Surface(500279) = {500278, 500394};
      Plane Surface(500281) = {500280};
      Plane Surface(500283) = {500282};
      Plane Surface(500285) = {500284};
      Plane Surface(500287) = {500286};
      Plane Surface(500289) = {500288};
      Plane Surface(500291) = {500290,500006,500234};
      Plane Surface(500295) = {500294};
      Plane Surface(500297) = {500296};
      Plane Surface(500299) = {500298};
      Plane Surface(500309) = {500308};
      Plane Surface(500311) = {500310};
      Plane Surface(500313) = {500312};
      Plane Surface(500315) = {500314};
      Plane Surface(500317) = {500316};
      Plane Surface(500319) = {500318};
      Plane Surface(500321) = {500320};
      Plane Surface(500323) = {500322};
      Plane Surface(500333) = {500332};
      Plane Surface(500335) = {500334};
      Plane Surface(500337) = {500336};
      Plane Surface(500343) = {500342,500356};
      Plane Surface(500345) = {500344};
      Plane Surface(500353) = {500352};
      Plane Surface(500359) = {500358,500372};
      Plane Surface(500361) = {500360};
      Plane Surface(500369) = {500368};
      Plane Surface(500375) = {500374};
      Plane Surface(500377) = {500376};
      Plane Surface(500379) = {500378};
      Plane Surface(500389) = {500388};
      Plane Surface(500391) = {500390};
      Plane Surface(500393) = {500392};
      Plane Surface(500399) = {500398};
      Plane Surface(500401) = {500400};
      Plane Surface(500407) = {500406};
      Plane Surface(500409) = {500408};
      Plane Surface(500415) = {500414};
      Plane Surface(500417) = {500416};
      Plane Surface(500421) = {500422};
      Plane Surface(500423) = {500424};
      Plane Surface(500425) = {500426};
      Plane Surface(500427) = {500430};
      Plane Surface(500429) = {500432};
      Plane Surface(500431) = {500434};
      Plane Surface(500433) = {500436};
      Plane Surface(500435) = {500438};
      Plane Surface(500437) = {500440};


      Field[1] = Attractor;
      Field[1].NNodesByEdge = 100;
      Field[1].EdgesList = {128,136,158,150,162,121,152};

      Field[2] = Attractor;
      Field[2].NNodesByEdge = 100;
      //Field[2].EdgesList = {294,296,289,277,284,273,282, 286};
      Field[2].EdgesList = {294,296,284,273,282, 286, 307,308,316,317,309,318};

      Field[6] = Attractor;
      Field[6].NNodesByEdge = 100;
      Field[6].EdgesList = {216,217,224,228};

      Field[3] = Threshold;
      Field[3].IField = 1;
      Field[3].LcMin = lc / 30;
      Field[3].LcMax = lc;
      Field[3].DistMin = 0.15;
      Field[3].DistMax = 10;

      Field[4] = Threshold;
      Field[4].IField = 2;
      Field[4].LcMin = lc / 20;
      Field[4].LcMax = lc;
      Field[4].DistMin = 0.15;
      Field[4].DistMax = 6;

      Field[7] = Threshold;
      Field[7].IField = 2;
      Field[7].LcMin = lc / 20;
      Field[7].LcMax = lc;
      Field[7].DistMin = 0.15;
      Field[7].DistMax = 6;

      Field[5] = Min;
      Field[5].FieldsList = {3,4,7};
      Background Field = 5;
      Mesh.Algorithm = 2;
    """

    cube_geometry = (
        "lc = " + str(h) + ";\n" + cube_stub)

    return __generate_grid_from_geo_string(cube_geometry)


def reentrant_cube_exterior(h=0.1, refinement_factor=0.3):
    reentrant_cube_stub = """
    Point(1) = {0, 0, 0, h};
    Point(2) = {1, 0, 0, h};
    Point(3) = {1, 1, 0, h};
    Point(4) = {0, 1, 0, h};
    Point(5) = {0, 0, 1, h};
    Point(6) = {1, 0, 1, h};
    Point(7) = {1, 1, 1, h};
    Point(8) = {0, 1, 1, h};
    Point(9) = {.5, .5, .5, h};
    Point(10) = {0, 0, .5, h};
    Point(11) = {0.5, 0, .5, h};
    Point(12) = {0, .5, .5, h};
    Point(13) = {0.5, 0, 1, h};
    Point(14) = {0.5, 0.5, 1, h};
    Point(15) = {0, 0.5, 1, h};

    Line(1) = {1, 2};
    Line(2) = {2, 6};
    Line(3) = {6, 13};
    Line(4) = {13, 11};
    Line(5) = {11, 10};
    Line(6) = {10, 1};

    Line(7) = {2, 3};
    Line(8) = {3, 7};
    Line(9) = {7, 6};
    
    Line(10) = {3, 4};
    Line(11) = {4, 8};
    Line(12) = {8, 7};

    Line(13) = {4, 1};
    Line(14) = {10, 12};
    Line(15) = {12, 15};
    Line(16) = {15, 8};

    Line(17) = {12, 9};
    Line(18) = {9, 14};
    Line(19) = {14, 15};

    Line(20) = {11, 9};
    Line(21) = {13, 14};

    Line Loop(1) = {1, 2, 3, 4, 5, 6};
    Line Loop(2) = {-2, 7, 8, 9};
    Line Loop(3) = {-8, 10, 11, 12};
    Line Loop(4) = {13, -6, 14, 15, 16, -11};
    Line Loop(5) = {17, 18, 19, -15};
    Line Loop(6) = {-14, -5, 20, -17};
    Line Loop(7) = {-4, 21, -18, -20};
    Line Loop(8) = {-13, -10, -7, -1};
    Line Loop(9) = {-21, -3, -9, -12, -16, -19};

    Plane Surface(1) = {1};
    Plane Surface(2) = {2};
    Plane Surface(3) = {3};
    Plane Surface(4) = {4};
    Plane Surface(5) = {5};
    Plane Surface(6) = {6};
    Plane Surface(7) = {7};
    Plane Surface(8) = {8};
    Plane Surface(9) = {9};

    Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    Volume(1) = {1};

    Field[1] = Attractor;
    Field[1].NodesList = {1,2,3,4,10,13,8,7,15,6};

    Field[2] = Threshold;
    Field[2].IField = 1;
    Field[2].LcMin = r;
    Field[2].LcMax = 3.5*h;
    Field[2].DistMin = 0.20;
    Field[2].DistMax = 0.20;

    Field[3] = Attractor;
    Field[3].NNodesByEdge = 100;
    //Field[3].EdgesList = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,21};
    Field[3].EdgesList = {1,2,3,4};
    Field[4] = Threshold;
    Field[4].IField = 3;
    Field[4].LcMin = 2*r;
    Field[4].LcMax = 4*h;
    Field[4].DistMin = 0.01;
    Field[4].DistMax = 0.01;

    Field[7] = Min;
    Field[7].FieldsList = {2};

    Background Field = 7;

    Mesh.Algorithm = 6;
    """
    reentrant_cube_geometry = (
        "h = " + str(h) + ";\n" +
        "r = h * " + str(refinement_factor) + ";\n" + reentrant_cube_stub)
    return __generate_grid_from_geo_string(reentrant_cube_geometry)


