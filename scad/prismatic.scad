tolerance = 0.015;
hole_radius = 0.045;
grid_hole_radius = 0.03;
outer_radius = 0.3;
thickness = 0.05;
hole_attach_height = 0.02;
compressed_height = 0.7;
extended_height = 1;
next_hole_radius = 0.045;
next_hole_attach_height = 0.02;
next_inner = 0.2;
hole_twist = 0;

//other params
eps = tolerance/100;
rail_attach_height = grid_hole_radius + hole_attach_height;
inner_radius = outer_radius - thickness;
rail_attach_width = rail_attach_height;
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

module bottom_joint() {
    out_hole_angle = 60;

    //grid of holes
    //hole_grid(hole_radius*2 + hole_attach_height * 2, grid_hole_radius);
   
    difference() {
        //outer cylinder
        cylinder(h=compressed_height, r= outer_radius, $fn = fnh);
        
        //shell
        translate([0,0,-eps/2])
        cylinder(h=compressed_height + eps, r= inner_radius, $fn = fnh);
        
        //screw holes
        translate([0,0,hole_attach_height + hole_radius]) 
        rotate([0,0,hole_twist]){
            rotate([90,0,0])
            cylinder(h=outer_radius*2 + eps,r=hole_radius,$fn=fnh,center=true);
            rotate([0,90,0])
            cylinder(h=outer_radius*2 + eps,r=hole_radius,$fn=fnh,center=true);
            
            
        }
        
        //space between guide holes and grid
        guide_hole_space = grid_hole_radius;
        //guide holes
        translate([-grid_hole_radius*cos(out_hole_angle) + inner_radius,0,hole_attach_height*2 + hole_radius*2 + grid_hole_radius + grid_hole_radius*sin(out_hole_angle) + guide_hole_space])
                rotate([0,out_hole_angle,0])
                cylinder(h=thickness/sin(out_hole_angle) + grid_hole_radius*2/tan(out_hole_angle), r=grid_hole_radius, $fn = fnh);
        
        translate([grid_hole_radius*cos(out_hole_angle) - inner_radius,0,hole_attach_height*2 + hole_radius*2 + grid_hole_radius + grid_hole_radius*sin(out_hole_angle) + guide_hole_space])
                rotate([0,-out_hole_angle,0])
                cylinder(h=thickness/sin(out_hole_angle) + grid_hole_radius*2/tan(out_hole_angle), r=grid_hole_radius, $fn = fnh);
        
        //for future reference of minimum joint size
        guide_height = hole_attach_height*2 + hole_radius*2 + grid_hole_radius + grid_hole_radius*2/sin(out_hole_angle) + thickness/tan(out_hole_angle) + guide_hole_space;
        
        //track
        rail_height = rail_attach_height + extended_height - compressed_height;
        translate([inner_radius,0,compressed_height - rail_height/2 - rail_attach_height])
            cube([thickness,rail_attach_width,rail_height],center = true);
        
        translate([-inner_radius,0,compressed_height - rail_height/2 - rail_attach_height])
            cube([thickness,rail_attach_width,rail_height],center = true);
        
        translate([0,0,compressed_height - rail_attach_height*2])
        cylinder(h=rail_attach_height,r=inner_radius + thickness/2,$fn=fn);
        
        //insert to track
        translate([0,0,compressed_height - rail_attach_height])
        cube([rail_attach_width, inner_radius*2 + thickness, rail_attach_height*2 + eps],center=true);
        
        //groove for string to slide through
        translate([0,0,compressed_height])
        rotate([0,90,0])
        cylinder(h=outer_radius*2 + eps, r=grid_hole_radius,center=true, $fn=fnh);
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
    difference() {
        union() {
            //loop for strings
            translate([0, inner_radius - tolerance - grid_hole_radius - hole_attach_height,0])
            cylinder(h=rail_attach_height - tolerance,r=grid_hole_radius + hole_attach_height,$fn=fnh);
            
            translate([0, -inner_radius + tolerance + grid_hole_radius + hole_attach_height,0])
            cylinder(h=rail_attach_height - tolerance,r=grid_hole_radius + hole_attach_height,$fn=fnh);
        }
        
        translate([0,0,-eps/2]) {
            //hole for string loops
        translate([0, inner_radius - tolerance - grid_hole_radius - hole_attach_height,0])
            cylinder(h=rail_attach_height-tolerance + eps,r=grid_hole_radius,$fn=fnh);
            
        translate([0, -inner_radius + tolerance + grid_hole_radius + hole_attach_height,0])
        cylinder(h=rail_attach_height - tolerance + eps,r=grid_hole_radius,$fn=fnh);
        }
    }
    
    
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
        translate([0,0,next_hole_attach_height + next_hole_radius + total_height]) {
            angled_cylinder(0, (next_inner + thickness + eps) * 2,  next_hole_radius, 90);
            angled_cylinder(90, (next_inner + thickness + eps) * 2, next_hole_radius, 0);
        }
    }
 
}

bottom_joint();

translate([outer_radius + 0.5, 0, 0])
top_joint();