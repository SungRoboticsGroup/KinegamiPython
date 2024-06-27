tolerance=0.075;
hole_radius=0.05;
grid_hole_radius=0.05;
outer_radius=1;
thickness=0.2;
hole_attach_height=0.05;
compressed_height=3.3134;
extended_height=5.8654;
next_hole_radius=0.05;
next_hole_attach_height=0.05;
next_inner=0.8;
hole_twist=0.0017;
//other params
eps = tolerance/100;
rail_attach_height = grid_hole_radius + hole_attach_height;
inner_radius = outer_radius - thickness;
rail_attach_width = rail_attach_height;
fn = 12;
fnh = 6;

module angled_cylinder(angle_deg_y, _h, _r, angle_deg_x = 0) {
    rotate([angle_deg_x, angle_deg_y, 0]) { // Rotate the cylinder 
            cylinder(h = _h, r = _r, center = true, $fn = fnh);
        }
}

module top_joint() {
    total_height = extended_height - compressed_height + rail_attach_height;
    difference() {
        union() {
            //outer cylinder
            cylinder(h=extended_height - compressed_height + rail_attach_height, r = inner_radius - tolerance, $fn = fn);
            
            //rail attachments
            translate([-(rail_attach_width - tolerance)/2,inner_radius - tolerance - thickness/2, 0])
            cube([rail_attach_width - tolerance, thickness, rail_attach_height - tolerance]);
            
            translate([-(rail_attach_width - tolerance)/2,-inner_radius + tolerance - thickness/2, 0])
            cube([rail_attach_width - tolerance, thickness, rail_attach_height - tolerance]);
            
        }
        //shell out inner cylinder
        translate([0,0,-eps/2])
        cylinder(h=extended_height - compressed_height + rail_attach_height + eps, r= inner_radius - tolerance - thickness, $fn = fn);
        
        translate([0,0,rail_attach_height + (extended_height - compressed_height)/2])
        cube([grid_hole_radius*2, inner_radius*2 + eps, extended_height - compressed_height + eps],center=true);
    }
    
    //to tie strings
//    difference() {
//        union() {
//            //loop for strings
//            translate([0, inner_radius - tolerance - grid_hole_radius - hole_attach_height,0])
//            cylinder(h=rail_attach_height - tolerance,r=grid_hole_radius + hole_attach_height,$fn=fnh);
//            
//            translate([0, -inner_radius + tolerance + grid_hole_radius + hole_attach_height,0])
//            cylinder(h=rail_attach_height - tolerance,r=grid_hole_radius + hole_attach_height,$fn=fnh);
//        }
//        
//        translate([0,0,-eps/2]) {
//            //hole for string loops
//        translate([0, inner_radius - tolerance - grid_hole_radius - hole_attach_height,0])
//            cylinder(h=rail_attach_height-tolerance + eps,r=grid_hole_radius,$fn=fnh);
//            
//        translate([0, -inner_radius + tolerance + grid_hole_radius + hole_attach_height,0])
//        cylinder(h=rail_attach_height - tolerance + eps,r=grid_hole_radius,$fn=fnh);
//        }
//    }
    
    
    //attachment segment
    attach_thickness = hole_attach_height;
    difference() {
        union() {
            //attachment for smaller cylinder
            translate([0,0,-attach_thickness*2 + total_height])
            cylinder(h=attach_thickness*2, r=inner_radius - tolerance, $fn=fn);
            translate([0,0,-attach_thickness*2 + total_height])
            cylinder(h=attach_thickness*2, r1=inner_radius - tolerance, r2=next_inner - tolerance - thickness, $fn=fn);
            
            translate([0,0,total_height])
            //small cylinder
            cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2, r = next_inner - tolerance, $fn=fn);
        }
        
        //shelling
        translate([0,0, -attach_thickness*2 - eps/2 + total_height])
        cylinder(h=next_hole_radius*2 + next_hole_attach_height * 2 + attach_thickness*2 + eps, r = min(inner_radius, next_inner - tolerance - thickness), $fn=fn);
        
        //holes for screws
//        translate([0,0,next_hole_attach_height + next_hole_radius + total_height]) {
//            angled_cylinder(0, (next_inner + thickness + eps) * 2,  next_hole_radius, 90);
//            angled_cylinder(90, (next_inner + thickness + eps) * 2, next_hole_radius, 0);
//        }
    }
 
}

translate([0,0,-extended_height+compressed_height-rail_attach_height])
top_joint();