tolerance=0.01;
hole_radius=0.05;
grid_hole_radius=0.1;
outer_radius=1;
thickness=0.2;
hole_attach_height=0.05;
attach_thickness=0.2;

//parameters that don't have to be defined explicitly
eps = tolerance/100;
inner_radius = outer_radius - thickness;
fn = outer_radius * 100;
fnh = hole_radius * 1000;

module angled_cylinder(angle_deg_y, _h, _r, angle_deg_x = 0) {
    rotate([angle_deg_x, angle_deg_y, 0]) { // Rotate the cylinder 
            cylinder(h = _h, r = _r, center = true, $fn = fnh);
        }
}

module hole(_h) {
    cylinder(h = _h, r = grid_hole_radius, $fn = fnh);
}

//makes circle of holes
module circular_grid(_h, _r, _n) {
    for (i = [0:_n-1]) {
        angle = i * 360 / _n;
        x_pos = _r * cos(angle);
        y_pos = _r * sin(angle);
        translate([x_pos, y_pos, -eps]) {
            hole(_h);
        }
    }
}

//makes grid of holes
module hole_grid(hole_height, grid_h) {
    translate([0,0,hole_height]) {
        difference() {
            //cylinder
            cylinder(h = grid_h, r = (inner_radius + outer_radius)/2, 
                center=false, $fn = fn);
            //holes
            for (i = [0 : grid_hole_radius * 2 + hole_attach_height : inner_radius - grid_hole_radius * 2]) {
                circular_grid(grid_h + eps*2, i, ceil(6.28*i/(grid_hole_radius*4)));
            }
            
        } 
    }
}

//makes curve of CSC
module curve(_radius, _twist, _turning_angle, _turning_radius) {
    rotate([90,0,_twist])
    translate([-_turning_radius, 0, 0])
    rotate_extrude(angle=_turning_angle, convexity=100, $fn = fn)
    translate([_turning_radius, 0]) circle(_radius);
}

module CSC(_r, _t1, _ta1, _tr1, _h, _t2, _ta2, _tr2, is_inner) {
    
    if (_ta1 > 90) {
        CSC(_r, _t1, 90, _tr1, 0, 0, 0, 0, is_inner);
        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1])
        rotate([0,-90,_t1])
        CSC(_r, 0, _ta1 - 90, _tr1, _h, _t2, _ta2, _tr2, is_inner);
    } else if (_ta2 > 90) {
        CSC(_r, _t1, _ta1, _tr1, _h, _t2, 90, _tr2, is_inner);
        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        translate([pos2*cos(_t2),pos2*sin(_t2),_tr2])
        rotate([0,-90,_t2])
        CSC(_r, 0, 0, 0, 0, 0, _ta2 - 90, _tr2, is_inner);
    } else {
        //intermediate positions
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2*sin(_ta2),2)));
        
        union() {
            //bottom shell threshold
            if (is_inner) {
                cylinder(h=eps, r=_r, center=true);
            }
            
            //first curve of CSC
            curve(_r, _t1, _ta1, _tr1);
            
            //straight link of CSC
            translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
            rotate([0,-_ta1,_t1])
            if (is_inner) {
                //add threshold for shelling
                translate([0,0,-eps /2])
                cylinder(h=_h + eps,r=_r, center = false, $fn= fn);
            } else {
                cylinder(h=_h,r=_r, center = false, $fn= fn);
            }
            
            //second curve of CSC
            translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
            rotate([0,-_ta1,_t1])
            translate([0, 0, _h])
            curve(_r, _t2, _ta2, _tr2);
            
            //shelling thresholds
            if (is_inner) {
                translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
                rotate([0,-_ta1,_t1])
                translate([0, 0, _h])
                translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
                rotate([0,-_ta2,_t2]) {
                    cylinder(h=eps, r=_r, center=true,$fn=fn);
                }
            
            }
        }
    }
    
}

//allow modules to attach
module CSC_ATTACH(_r, _t1, _ta1, _tr1, _h, _t2, _ta2, _tr2, next_hole_radius, next_hole_attach_height, next_inner) {
    
    if (_ta1 > 90) {        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1])
        rotate([0,-90,_t1])
        CSC_ATTACH(_r, 0, _ta1 - 90, _tr1, _h, _t2, _ta2, _tr2, next_hole_radius, next_hole_attach_height, next_inner);
    } else if (_ta2 > 90) {        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        translate([pos2*cos(_t2),pos2*sin(_t2),_tr2])
        rotate([0,-90,_t2])
        CSC_ATTACH(_r, 0, 0, 0, 0, 0, _ta2 - 90, _tr2, next_hole_radius, next_hole_attach_height, next_inner);
    } else {
        //intermediate positions
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2*sin(_ta2),2)));
        
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
        rotate([0,-_ta2,_t2]) {
            difference() {
                union() {
                    //attachment for smaller cylinder
                    translate([0,0,-attach_thickness])
                    cylinder(h=attach_thickness, r=_r - thickness/2, $fn=fn);
                    translate([0,0,-attach_thickness*2])
                    cylinder(h=attach_thickness*2, r1=_r - thickness/2, r2=next_inner - tolerance - thickness, $fn=fn);
                    
                    //small cylinder
                    cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2, r = next_inner - tolerance, $fn=fn);
                }
                
                //shelling
                translate([0,0, -attach_thickness*2 - eps/2])
                cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2 + attach_thickness*2 + eps, r = min(inner_radius, next_inner - tolerance - thickness), $fn=fn);
                
                //holes for screws
                translate([0,0,next_hole_attach_height + next_hole_radius]) {
                    angled_cylinder(0, (next_inner + thickness + eps) * 2,  hole_radius, 90);
                    angled_cylinder(90, (next_inner + thickness + eps) * 2, hole_radius, 0);
                }
            }
        }
    }
}

//remove overlap to allow all branches to attach to each other
module CSC_REMOVE_OVERLAP(_r, _t1, _ta1, _tr1, _h, _t2, _ta2, _tr2, next_hole_radius, next_hole_attach_height, next_inner) {
    if (_ta1 > 90) {        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1])
        rotate([0,-90,_t1])
        CSC_REMOVE_OVERLAP(_r, 0, _ta1 - 90, _tr1, _h, _t2, _ta2, _tr2, next_hole_radius, next_hole_attach_height, next_inner);
    } else if (_ta2 > 90) {        
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2,2)));
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        translate([pos2*cos(_t2),pos2*sin(_t2),_tr2])
        rotate([0,-90,_t2])
        CSC_REMOVE_OVERLAP(_r, 0, 0, 0, 0, 0, _ta2 - 90, _tr2, next_hole_radius, next_hole_attach_height, next_inner);
    } else {
        pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
        pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2*sin(_ta2),2)));
        
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
        rotate([0,-_ta2,_t2]) {
//            cylinder(h=999999, r=_r + tolerance, $fn=fn);
        }
    }
}

//c is _t1, _ta1, _tr1, _h, _t2, _ta2, _tr2, next hole radius, next hole height, next inner radius
module branch(branch_list, _outer, _inner) {
    union() {
        //make grid of holes
        //hole_grid(hole_radius*2 + hole_attach_height * 2, grid_hole_radius);
        
        //make CSC shells for each CSC in branch_list
        difference() {
            for (c = branch_list) {
                CSC(_outer,c[0],c[1],c[2],c[3],c[4],c[5],c[6], false);
            }
            
            for (c = branch_list) {
                CSC(_inner,c[0],c[1],c[2],c[3],c[4],c[5],c[6], true);
                CSC_REMOVE_OVERLAP(_outer,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7], c[8], c[9]);
            }
            
            //make screw holes at the bottom
            translate([0,0,hole_radius + hole_attach_height]) {
                angled_cylinder(0, (outer_radius*1.5 + eps) * 2,  hole_radius, 90);
                angled_cylinder(90, (outer_radius*1.5 + eps) * 2, hole_radius, 0);
            }    
        }
        for (c=branch_list) {     
            CSC_ATTACH(_outer,c[0],c[1],c[2],c[3],c[4],c[5],c[6], c[7], c[8], c[9]);
        }
    }
}

branch([
[ 0.0, 42.750838323, 1.0*outer_radius, 1.545671284*outer_radius, -0.0, 2.249161677, 1.0*outer_radius, 0.05, 0.05, 0.8],
[ 180.0, 0.0, 1.0*outer_radius, 6.3625*outer_radius, 180.0, 0.0, 1.0*outer_radius, 0.05, 0.05, 0.8],
[ 180.0, 63.377666873, 1.0*outer_radius, 1.282149075*outer_radius, -180.0, 63.377666873, 1.0*outer_radius, 0.05, 0.05, 0.8],
],outer_radius,inner_radius);