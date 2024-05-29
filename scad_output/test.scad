eps=0.0003;
hole_radius=0.045;
grid_hole_radius=0.03;
outer_radius=0.3;
inner_radius=0.25;
fn = outer_radius * 100;
fnh = hole_radius * 1000;

module angled_cylinder(angle_deg_y, _h, _r, angle_deg_x = 0) {
    rotate([angle_deg_x, angle_deg_y, 0]) { // Rotate the cylinder 
            cylinder(h = _h, r = _r, center = true, $fn = fnh);
        }
}

module hole(_h) {
    cylinder(h = _h, r = grid_hole_radius, $fn = fnh); // Adjust height (h) as needed
}

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

module hole_grid(hole_height, grid_h) {
    difference() {
        translate([0,0,hole_height - hole_radius])
        cylinder(h = grid_h, r = (inner_radius + outer_radius)/2, center=false, $fn = fn);
        translate([0,0,hole_height - hole_radius])
        for (i = [0 : hole_radius * 2 : inner_radius - grid_hole_radius * 2]) {
            circular_grid(grid_h + eps*2, i, ceil(6.28*i/(grid_hole_radius*4)));
        }
    }
}

module curve(_radius, _twist, _turning_angle, _turning_radius) {
    rotate([90,0,_twist])
    translate([-_turning_radius, 0, 0])
    rotate_extrude(angle=_turning_angle, convexity=100, $fn = fn)
       translate([_turning_radius, 0]) circle(_radius);
}

module CSC(_r, _t1, _ta1, _tr1, _h, _t2, _ta2, _tr2, is_inner) {
    pos1 = (-_tr1 + sqrt(_tr1*_tr1 - pow(_tr1*sin(_ta1),2)));
    pos2 = (-_tr2 + sqrt(_tr2*_tr2 - pow(_tr2*sin(_ta2),2)));
    union() {
        if (is_inner) {
            cylinder(h=eps, r=_r, center=true);
        }
        
        curve(_r, _t1, _ta1, _tr1);
        
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        if (is_inner) {
            translate([0,0,-eps /2])
            cylinder(h=_h + eps,r=_r, center = false, $fn= fn);
        } else {
            cylinder(h=_h,r=_r, center = false, $fn= fn);
        }
        
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
        rotate([0,-_ta1,_t1])
        translate([0, 0, _h])
        //cylinder(h=_h - 4,r=_r + 2,center = false);
        curve(_r, _t2, _ta2, _tr2);
        
        if (is_inner) {
            translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
            rotate([0,-_ta1,_t1])
            translate([0, 0, _h])
            translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
            rotate([0,-_ta2,_t2])
            cylinder(h=eps, r=_r, center=true,$fn=fn);
            
            translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
            rotate([0,-_ta1,_t1])
            translate([0, 0, _h])
            translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
            rotate([0,-_ta2,_t2])
            translate([0,0,-hole_radius*1.5])
        angled_cylinder(0, (outer_radius + eps) * 2.5,  hole_radius, 90);
        translate([pos1*cos(_t1),pos1*sin(_t1),_tr1*sin(_ta1)])
            rotate([0,-_ta1,_t1])
            translate([0, 0, _h])
            translate([pos2*cos(_t2),pos2*sin(_t2),_tr2*sin(_ta2)])
            rotate([0,-_ta2,_t2])
        translate([0,0,-hole_radius*1.5])
    angled_cylinder(90, (outer_radius + eps) * 2.5, hole_radius, 0);
        }
    }
}

//_t1, _ta1, _tr1, _h, _t2, _ta2, _tr2
module branch(branch_list, _outer, _inner) {
    union() {
        hole_grid(hole_radius*4, hole_radius*2);
        difference() {
        union() {
            for (c = branch_list) {
                CSC(_outer,c[0],c[1],c[2],c[3],c[4],c[5],c[6], false);
            }
         }
        union() {
            for (c = branch_list) {
                CSC(_inner,c[0],c[1],c[2],c[3],c[4],c[5],c[6], true);
            }
        }
        translate([0,0,hole_radius * 1.5])
        angled_cylinder(0, (outer_radius*1.5 + eps) * 2,  hole_radius, 90);
        translate([0,0,hole_radius * 1.5])
    angled_cylinder(90, (outer_radius*1.5 + eps) * 2, hole_radius, 0);
    }
    }
}

branch([[ 0.0, 39.101192967, 1.0*outer_radius, 2.310330337*outer_radius, -0.0, 5.898807033, 1.0*outer_radius],
[ 0.0, 0.0, 1.0*outer_radius, 7.292893219*outer_radius, 180.0, 0.0, 1.0*outer_radius],
[ 180.0, 41.615509189, 1.0*outer_radius, 2.627859538*outer_radius, 180.0, 41.615509189, 1.0*outer_radius],
],outer_radius,inner_radius);