#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform vec4 u_Color;
uniform float u_Anim;

in vec2 fs_Pos;
out vec4 out_Col;

vec3 light;
vec3 nor;
float worl;

vec3 anim_trans;
vec3 anim_angle;

float map_value;
bool water_hit;

vec3 rayCast(vec4 s) {
    // multiply by far clip plane
    float far_clip = 1000.0;
    float near_clip = 0.1;
    s *= far_clip;

    // multiply by inverse projection matrix
    mat4 proj;
    float fov = 45.0;
    float fov_rad = (fov * 3.14159) / 180.0;
    float S_f = tan(fov_rad / 2.0);
    float a = 1.0 / ((u_Dimensions.x / u_Dimensions.y) * S_f);
    float b = 1.0 / S_f;
    float P = far_clip / (far_clip - near_clip);
    float Q = (-1.f * far_clip * near_clip) / (far_clip - near_clip);
    proj[0][0] = a;
    proj[1][1] = b;
    proj[2][2] = P;
    proj[3][2] = Q;
    proj[2][3] = 1.0;
    proj[3][3] = 0.0;
    s = inverse(proj) * s;

    // multiply by inverse of view matrix
    mat4 view;
    mat4 orient;
    mat4 transl;
    vec3 forw_axis = u_Ref - u_Eye;
    vec3 right_axis = cross(u_Up, forw_axis);
    vec3 up_axis = cross(right_axis, forw_axis);

    // special case where forward is world up
    if (forw_axis == u_Up) {
      right_axis = vec3(1.0, 0.0, 0.0);
      up_axis = vec3(0.0, 0.0, -1.0);
    }

    // set up orientation and translation matrices
    for (int i = 0; i < 3; i++) {
      orient[0][i] = right_axis[i];
      orient[1][i] = up_axis[i];
      orient[2][i] = forw_axis[i];
      transl[3][i] = u_Eye[i] * -1.0;
    }
    view = orient * transl;
    s = inverse(view) * s;

    // set up ray
    vec3 origin = u_Eye;
    vec3 dir = normalize(vec3(s) - u_Eye);

    // set light vector for shading
    light = vec3((inverse(view) * vec4(0.0, 0.0, 0.0, 1.0)) - vec4(fs_Pos, 1.0, 1.0));

    return dir;
}

float random1(float t) {
    return 2.0 * fract(sin(t * 489.12342) * 348921.32457) - 1.0;
}

vec2 random2(vec2 p, vec2 seed) {
  return fract(sin(vec2(dot(p + seed, vec2(311.7, 127.1)), dot(p + seed, vec2(269.5, 183.3)))) * 85734.3545);
}

vec3 remapColor(vec3 c, float offset, float seed) {
    float ro = random1(seed) * offset;
    float go = random1(seed + 31242.134) * offset;
    float bo = random1(seed + 73576.347) * offset;
    
    return clamp(vec3(c.r + ro, c.g + go, c.b + bo), 0.0, 1.0);
}

float WorleyNoise(vec2 uv) {
    // Tile the space
    vec2 uvInt = floor(uv);
    vec2 uvFract = fract(uv);

    float minDist = 1.0; // Minimum distance initialized to max.

    // Search all neighboring cells and this cell for their point
    for(int y = -1; y <= 1; y++) {
        for(int x = -1; x <= 1; x++) {
            vec2 neighbor = vec2(float(x), float(y));

            // Random point inside current neighboring cell
            vec2 point = random2(uvInt + neighbor, vec2(10.0));

            // Compute the distance b/t the point and the fragment
            // Store the min dist thus far
            vec2 diff = neighbor + point - uvFract;
            float dist = length(diff);
            minDist = min(minDist, dist);
        }
    }
    return minDist;
}

float sdf_sphere(vec3 p, float rad) {
  float sphere = length(p) - rad;
  return sphere;
}

float sdf_box(vec3 p, vec3 b) {
  vec3 d = abs(p) - b;
  return length(max(d,0.0));
}

float sdf_torus(vec3 p, vec2 t) {
  vec2 q = vec2(length(p.xz) - t.x, p.y);
  return length(q) - t.y;
}

float sdf_cylin(vec3 p, vec2 h) {
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float dot2(vec2 v) {
  return dot(v, v); 
}

float sdf_capcone(vec3 p, float h, float r1, float r2) {
  vec2 q = vec2(length(p.xz), p.y);
  vec2 k1 = vec2(r2, h);
  vec2 k2 = vec2(r2-r1, 2.0 * h);
  vec2 ca = vec2(q.x - min(q.x,(q.y < 0.0) ? r1:r2), abs(q.y) - h);
  vec2 cb = q - k1 + k2 * clamp(dot(k1 - q, k2) / dot2(k2), 0.0, 1.0);
  float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
  return s * sqrt(min(dot2(ca), dot2(cb)));
}

float sdf_vertcap(vec3 p, float h, float r) {
  p.y -= clamp(p.y, 0.0, h);
  return length(p) - r;
}

float elongcap_op(vec3 p, vec3 h, float ht, float r) {
  vec3 q = abs(p) - h;
  return sdf_vertcap(max(q, 0.0), ht, r) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float elongtor_op(vec3 p, vec3 h, vec2 t) {
  vec3 q = abs(p) - h;
  return sdf_torus(max(q, 0.0), t) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float union_op(float d1, float d2) {
  return min(d1, d2);
}

float sm_union_op(float d1, float d2, float k) {
  float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
  return mix(d2, d1, h) - k * h * (1.0 - h); 
}

float sub_op(float d1, float d2) {
  return max(-1.0 * d1, d2);
}

float sect_op(float d1, float d2) {
  return max(d1, d2);
}

mat4 rot_matrix(vec3 r) {
  // convert angle to radians
  r = (r * 3.14159) / 180.0;

  mat4 m;
  m[0][0] = 1.0;
  m[1][1] = 1.0;
  m[2][2] = 1.0;
  m[3][3] = 1.0;
  mat4 rot_x = m;
  rot_x[1][1] = cos(r.x);
  rot_x[1][2] = sin(r.x);
  rot_x[2][1] = -1.0 * sin(r.x);
  rot_x[2][2] = cos(r.x);
  mat4 rot_y = m;
  rot_y[0][0] = cos(r.y);
  rot_y[2][2] = cos(r.y);
  rot_y[2][0] = sin(r.y);
  rot_y[0][2] = -1.0 * sin(r.y);
  mat4 rot_z = m;
  rot_z[0][0] = cos(r.z);
  rot_z[1][1] = cos(r.z);
  rot_z[1][0] = -1.0 * sin(r.z);
  rot_z[0][1] = sin(r.z);

  return rot_x * rot_y * rot_z;
}

vec3 trans_pt(vec3 p, vec3 t) {
  return p + t;
}

vec3 rot_op(vec3 r, vec3 p) {
  mat4 tmat = rot_matrix(-1.0 * r);
  vec3 p_rot = vec3(tmat * vec4(p, 1.0));
  return p_rot;
}

float pedestalSDF(vec3 p) {
  // BASE
  // translate to make pedestal layering
  vec3 t_cyl = vec3(0.0, 3.5, 0.0);
  vec3 t_b2 = vec3(0.0, -0.1, 0.0);
  vec3 t_b3 = vec3(0.0, -0.2, 0.0);
  vec3 p_b2 = trans_pt(p, t_b2);
  vec3 p_b3 = trans_pt(p, t_b3);

  // translate entire pedestal down
  vec3 p_b1 = trans_pt(p, t_cyl);
  p_b2 = trans_pt(p_b2, t_cyl);
  p_b3 = trans_pt(p_b3, t_cyl);

  // cylinder SDFs
  float base1 = sdf_cylin(p_b1, vec2(1.6, 0.07));
  float base2 = sdf_cylin(p_b2, vec2(1.4, 0.07));
  float base3 = sdf_cylin(p_b3, vec2(1.2, 0.07));

  // combine shapes
  float dist_base = union_op(base1, base2);
  dist_base = union_op(dist_base, base3);

  // MID 
  // translate to make layering
  vec3 t_m1 = vec3(0.0, -0.3, 0.0);
  vec3 t_m2 = vec3(0.0, -0.6, 0.0);
  vec3 t_m3 = vec3(0.0, -1.0, 0.0);
  vec3 p_m1 = trans_pt(p, t_m1);
  vec3 p_m2 = trans_pt(p, t_m2);
  vec3 p_m3 = trans_pt(p, t_m3);

  // translate down
  p_m1 = trans_pt(p_m1, t_cyl);
  p_m2 = trans_pt(p_m2, t_cyl);
  p_m3 = trans_pt(p_m3, t_cyl);
  
  // cap cone and cylinder SDFs
  float mid1 = sdf_capcone(p_m1, 0.3, 1.0, 0.8);
  float mid2 = sdf_cylin(p_m2, vec2(0.8, 0.3));
  float mid3 = sdf_capcone(p_m3, 0.1, 0.8, 1.0);

  float dist_mid = sm_union_op(mid1, mid2, 0.05);
  dist_mid = sm_union_op(dist_mid, mid3, 0.05);

  // combine sections
  float dist = union_op(dist_base, dist_mid);

  return dist;
}

float ped2SDF(vec3 p) {
  // TOP
  // translate to make pedestal layering
  vec3 t_cyl = vec3(0.0, 3.5, 0.0);
  vec3 t_t1 = vec3(0.0, -1.2, 0.0);
  vec3 t_t2 = vec3(0.0, -1.35, 0.0);
  vec3 t_t3 = vec3(0.0, -1.45, 0.0);
  vec3 p_t1 = trans_pt(p, t_t1);
  vec3 p_t2 = trans_pt(p, t_t2);
  vec3 p_t3 = trans_pt(p, t_t3);

  // translate entire pedestal down
  p_t1 = trans_pt(p_t1, t_cyl);
  p_t2 = trans_pt(p_t2, t_cyl);
  p_t3 = trans_pt(p_t3, t_cyl);

  // cap cone and cylinder SDFs
  float top1 = sdf_capcone(p_t1, 0.1, 0.8, 0.6);
  float top2 = sdf_capcone(p_t2, 0.08, 0.6, 0.7);
  float top3 = sdf_cylin(p_t3, vec2(0.5, 0.1));

  // combine shapes
  float dist_top = sm_union_op(top1, top2, 0.05);
  dist_top = sm_union_op(dist_top, top3, 0.05);
  return dist_top;
}

float treeSDF(vec3 p) {
  // TRUNKS
  // translate to place on pedestal
  vec3 t_cyl = vec3(0.0, 3.5, 0.0);
  vec3 t_tr1 = vec3(0.0, -5.0, 0.0);
  vec3 t_tr2a = vec3(2.8, 2.3, 0.0);
  vec3 t_tr2b = vec3(1.8, -1.68, 0.0);
  vec3 t_tr3a = vec3(-1.35, 1.6, 0.0);
  vec3 t_tr3b = vec3(-0.7, -2.88, 0.0);
  vec3 t_c1 = vec3(0.05, -7.5, 0.0);
  vec3 t_s2 = vec3(3.5, -0.2, 0.0);
  vec3 t_s3 = vec3(-0.1, 0.5, 0.0);
  vec3 t_s4 = vec3(-2.1, -2.0, 0.0);
  vec3 p_tr1 = trans_pt(p, t_tr1);
  vec3 p_tr2a = trans_pt(p, t_tr2a);
  vec3 p_tr2b = trans_pt(p, t_tr2b);
  vec3 p_tr3a = trans_pt(p, t_tr3a);
  vec3 p_tr3b = trans_pt(p, t_tr3b);
  vec3 p_c1 = trans_pt(p, t_c1);
  vec3 p_s2 = trans_pt(p, t_s2);
  vec3 p_s3 = trans_pt(p, t_s3);
  vec3 p_s4 = trans_pt(p, t_s4);

  // rotation of torus
  vec3 r_tr2a = vec3(90.0, 0.0, 0.0);
  p_tr2a = rot_op(r_tr2a, p_tr2a);
  p_tr2b = rot_op(r_tr2a, p_tr2b);
  p_tr3a = rot_op(r_tr2a, p_tr3a);
  p_tr3b = rot_op(r_tr2a, p_tr3b);

  // translate entire tree down
  p_tr1 = trans_pt(p_tr1, t_cyl);
  p_c1 = trans_pt(p_c1, t_cyl);

  // cylinder SDFs
  float trunk1 = sdf_cylin(p_tr1, vec2(0.25, 3.75));

  // torus SDFs
  vec3 e_tr2a = vec3(0.5, 0.0, 0.0);
  vec3 e_tr2b = vec3(0.5, 0.0, 0.3);
  vec3 e_tr3a = vec3(0.0, 0.0, 0.8);
  vec3 e_tr3b = vec3(0.0, 0.0, 1.8);
  float trunk2a = elongtor_op(p_tr2a, e_tr2a, vec2(2.0, 0.08));
  float trunk2b = elongtor_op(p_tr2b, e_tr2b, vec2(1.7, 0.08));
  float trunk3a = elongtor_op(p_tr3a, e_tr3a, vec2(1.0, 0.08));
  float trunk3b = elongtor_op(p_tr3b, e_tr3b, vec2(1.0, 0.08));

  // cut shapes SDFs
  vec3 e_c1 = vec3(0.0, 0.0, 1.0);
  float cut1 = elongcap_op(p_c1, e_c1, 2.5, 0.04);

  // sect shapes SDFs
  float sect1 = sdf_box(p, vec3(2.0, 2.0, 1.0));
  float sect2 = sdf_box(p_s2, vec3(1.5, 1.3, 1.0));
  float sect3 = sdf_box(p_s3, vec3(1.0, 1.5, 1.0));
  float sect4 = sdf_box(p_s4, vec3(1.0, 3.0, 1.0));

  // combine shapes
  //float dist_tree = union_op(trunk1a, trunk1b);
  //dist_tree = union_op(trunk1, dist_tree);
  //float dist_tree = union_op(trunk1, cut1);
  float dist_tree = sub_op(cut1, trunk1);
  //dist_tree = union_op(dist_tree, trunk2a);
  float dist_tr2a = sect_op(trunk2a, sect1);
  float dist_tr2b = sect_op(trunk2b, sect2);
  

  float dist_tr3a = sect_op(trunk3a, sect3);
  float dist_tr3b = sect_op(trunk3b, sect4);
  float dist_tr3 = union_op(dist_tr3a, dist_tr3b);
  dist_tree = union_op(dist_tree, dist_tr3);

  // BRANCHES
  vec3 t_b1 = vec3(-0.2, -3.2, 0.0);
  vec3 t_sb1 = vec3(-1.0, -2.0, 0.0);
  vec3 p_b1 = trans_pt(p, t_b1);
  vec3 p_sb1 = trans_pt(p, t_sb1);

  vec3 t_b2 = vec3(0.0, -3.5, 0.0);
  vec3 t_sb2 = vec3(1.0, -2.0, 0.0);
  vec3 p_b2 = trans_pt(p, t_b2);
  vec3 p_sb2 = trans_pt(p, t_sb2);
  vec3 t_b3 = vec3(4.0, -1.5, 0.0);
  vec3 t_sb3 = vec3(3.5, -1.0, 0.0);
  vec3 p_b3 = trans_pt(p, t_b3);
  vec3 p_sb3 = trans_pt(p, t_sb3);
  vec3 t_sb4 = vec3(4.5, -1.0, 0.0);
  vec3 t_sb4b = vec3(4.0, -0.4, 0.0);
  vec3 p_sb4 = trans_pt(p, t_sb4);
  vec3 p_sb4b = trans_pt(p, t_sb4b);
  vec3 t_b5 = vec3(3.5, -1.5, 0.0);
  vec3 t_sb5 = vec3(3.3, -1.3, 0.0);
  vec3 p_b5 = trans_pt(p, t_b5);
  vec3 p_sb5 = trans_pt(p, t_sb5);

  vec3 t_b6 = vec3(-1.7, -4.0, 0.0);
  vec3 p_b6 = trans_pt(p, t_b6);
  vec3 t_sb6 = vec3(-2.7, -4.0, 0.0);
  vec3 p_sb6 = trans_pt(p, t_sb6);

  vec3 t_c2 = vec3(4.0, -0.8, 0.0);
  vec3 p_c2 = trans_pt(p, t_c2);

  // rotation of torus
  p_b1 = rot_op(r_tr2a, p_b1);
  p_b2 = rot_op(r_tr2a, p_b2);
  p_b3 = rot_op(r_tr2a, p_b3);
  p_b5 = rot_op(r_tr2a, p_b5);
  p_b6 = rot_op(r_tr2a, p_b6);

  // torus SDFs
  vec3 e_b3 = vec3(0.0, 0.0, 0.3);
  vec3 e_b5 = vec3(0.0, 0.0, 0.1);
  vec3 e_b6 = vec3(0.0, 0.0, 1.5);
  float branch1 = elongtor_op(p_b1, e_tr3b, vec2(1.0, 0.16));
  float branch2 = elongtor_op(p_b2, e_tr3b, vec2(0.7, 0.16));
  float branch3 = elongtor_op(p_b3, e_b3, vec2(0.5, 0.05));
  float branch4 = elongtor_op(p_b3, e_b3, vec2(0.8, 0.05));
  float branch5 = elongtor_op(p_b5, e_b5, vec2(0.2, 0.05));
  float branch6 = elongtor_op(p_b6, e_b6, vec2(0.6, 0.07));

  // cut shapes SDFs
  float cut2 = sdf_box(p_c2, vec3(1.5, 0.7, 0.2));

  // sect shapes SDFs
  float sect_b1 = sdf_box(p_sb1, vec3(1.0, 2.5, 1.0));
  float sect_b2 = sdf_box(p_sb2, vec3(1.0, 2.5, 1.0));
  float sect_b3 = sdf_box(p_sb3, vec3(0.3, 0.6, 0.2));
  float sect_b4 = sdf_box(p_sb4, vec3(0.5, 0.8, 0.2));
  float sect_b4b = sdf_box(p_sb4b, vec3(0.2, 0.2, 0.2));
  float sect_b5 = sdf_box(p_sb5, vec3(0.2, 0.3, 0.2));
  float sect_b6 = sdf_box(p_sb6, vec3(1.0, 2.5, 1.0));

  // combine shapes
  float dist_tr2 = union_op(dist_tr2a, dist_tr2b);
  float dist_br = sect_op(sect_b1, branch1);
  float dist_br2 = sect_op(sect_b2, branch2);
  float dist_br3 = sect_op(sect_b3, branch3);
  float dist_br5 = sect_op(sect_b5, branch5);
  float sect_4 = union_op(sect_b4, sect_b4b);
  float dist_br4 = sect_op(sect_4, branch4);
  float dist_br6 = sect_op(sect_b6, branch6);

  dist_br3 = sect_op(cut2, dist_br3);
  dist_br4 = sect_op(cut2, dist_br4);
  dist_br5 = sect_op(cut2, dist_br5);

  dist_tree = sm_union_op(dist_tree, dist_tr2, 0.2);
  float dist = sm_union_op(dist_tree, dist_br, 0.1);
  dist = sm_union_op(dist_br2, dist, 0.03);
  dist = sm_union_op(dist_br3, dist, 0.03);
  dist = sm_union_op(dist_br4, dist, 0.03);
  dist = sm_union_op(dist_br5, dist, 0.03);
  dist = sm_union_op(dist_br6, dist, 0.03);

  return dist;
}

float waterSDF(vec3 p) {
  // translate
  vec3 t_w1 = vec3(4.2, -1.7, 0.0);
  vec3 p_w1 = trans_pt(p, t_w1);
  vec3 t_w2 = vec3(0.0, -3.5, 0.0);
  vec3 p_w2 = trans_pt(p, t_w2);
  vec3 t_w3 = vec3(0.0, 3.5, 0.0);
  vec3 p_w3 = trans_pt(p, t_w3);

  // box SDFs
  float water1 = sdf_box(p_w1, vec3(1.0, 1.0, 1.0));
  float water2 = sdf_box(p_w2, vec3(1.5, 1.5, 1.5));
  float water3 = sdf_box(p_w3, vec3(15.0, 1.0, 5.0));

  // combine shapes
  float dist = union_op(water1, water2);
  dist = union_op(water3, dist);

  return dist;
}

vec3 estNormalPed(vec3 p) {
  // find normal of pedestal points
  float eps = 0.001;
  vec3 nor_c = vec3(pedestalSDF(vec3(p.x + eps, p.y, p.z)) - pedestalSDF(vec3(p.x - eps, p.y, p.z)),
                  pedestalSDF(vec3(p.x, p.y + eps, p.z)) - pedestalSDF(vec3(p.x, p.y - eps, p.z)),
                  pedestalSDF(vec3(p.x, p.y, p.z + eps)) - pedestalSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_c);
}

vec3 estNormalPed2(vec3 p) {
  // find normal of pedestal 2 points
  float eps = 0.001;
  vec3 nor_c = vec3(ped2SDF(vec3(p.x + eps, p.y, p.z)) - ped2SDF(vec3(p.x - eps, p.y, p.z)),
                  ped2SDF(vec3(p.x, p.y + eps, p.z)) - ped2SDF(vec3(p.x, p.y - eps, p.z)),
                  ped2SDF(vec3(p.x, p.y, p.z + eps)) - ped2SDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_c);
}

vec3 estNormalTree(vec3 p) {
  // find normal of tree points
  float eps = 0.001;
  vec3 nor_c = vec3(treeSDF(vec3(p.x + eps, p.y, p.z)) - treeSDF(vec3(p.x - eps, p.y, p.z)),
                  treeSDF(vec3(p.x, p.y + eps, p.z)) - treeSDF(vec3(p.x, p.y - eps, p.z)),
                  treeSDF(vec3(p.x, p.y, p.z + eps)) - treeSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_c);
}

vec3 estNormalWater(vec3 p) {
  // find normal of water points
  float eps = 0.001;
  vec3 nor_c = vec3(waterSDF(vec3(p.x + eps, p.y, p.z)) - waterSDF(vec3(p.x - eps, p.y, p.z)),
                  waterSDF(vec3(p.x, p.y + eps, p.z)) - waterSDF(vec3(p.x, p.y - eps, p.z)),
                  waterSDF(vec3(p.x, p.y, p.z + eps)) - waterSDF(vec3(p.x, p.y, p.z - eps)));
  return normalize(nor_c);
}

vec2 rayMarch(vec3 eye, vec3 dir) { 
  // rayMarch returns (t, object id)
  float t = 0.01;
  int max_steps = 1000;
  vec3 p = eye + t * dir;
  for (int i = 0; i < max_steps; i++) {
    p = eye + t * dir;

    float dist = pedestalSDF(p);
    float dist2 = ped2SDF(p);
    float dist3 = treeSDF(p);
    //float dist4 = waterSDF(p);

    nor = estNormalPed(p);

    if (dist < 0.00001) {
      // at pedestal surface
      map_value = dist;
      nor = estNormalPed(p);
      return vec2(t, 1.0);  
      // move along ray
      t += dist;
    }
    else if (dist2 < 0.00001) {
      // at pedestal 2 surface
      map_value = dist2;
      nor = estNormalPed2(p);
      return vec2(t, 2.0);
      // move along ray
      t += dist2;
    }
    else if (dist3 < 0.00001) {
      // at tree surface
      map_value = dist3;
      nor = estNormalTree(p);
      return vec2(t, 3.0);
      // move along ray
      t += dist3;
    }
    else {
      // increment by smallest distance
      float dist_min = min(dist, dist2);
      dist_min = min(dist_min, dist3);
      t += dist_min;
      map_value = dist_min;
    }
  
    if (t >= 1000.0) {
      // end
      return vec2(t, 0.0);
    }
  }
  t = 1000.0;
  return vec2(t, 0.0);
}

vec2 rayMarchWater(vec3 eye, vec3 dir) { 
  // rayMarch returns (t, object id)
  float t = 0.01;
  int max_steps = 1000;
  vec3 p = eye + t * dir;
  for (int i = 0; i < max_steps; i++) {
    p = eye + t * dir;

    float dist = waterSDF(p);

    nor = estNormalWater(p);

    if (dist < 0.00001) {
      // at water surface
      map_value = dist;
      nor = estNormalWater(p);
      return vec2(t, 1.0);  
      // move along ray
      t += dist;
    }
    else {
      // increment by smallest distance
      t += dist;
    }
  
    if (t >= 1000.0) {
      // end
      return vec2(t, 0.0);
    }
  }
  t = 1000.0;
  return vec2(t, 0.0);
}

bool rayBoxIntersection(vec3 origin, vec3 dir, vec3 min, vec3 max) {
  // check intersection of ray with cube for bounding box purposes
  float near = -1.0 * (1.0 / 0.0);
  float far = (1.0 / 0.0);
  float t0;
  float t1;
  for (int i = 0; i < 3; i++) {
    if (dir[i] == 0.0) {
      if (origin[i] < min[i] || origin[i] > max[i]) {
        return false;
      }
    }
    t0 = (min[i] - origin[i]) / dir[i];
    t1 = (max[i] - origin[i]) / dir[i];
    if (t0 > t1) {
      float temp = t0;
      t0 = t1;
      t1 = temp;
    }
    if (t0 > near) {
      near = t0;
    }
    if (t1 < far) {
      far = t1;
    }
  }
  if (near > far) {
    return false;
  }
  else {
    return true;
  }
}

float softShadow(vec3 dir, vec3 origin, float min_t, float k) {
  float res = 1.0;
  float t = min_t;
  for(int i = 0; i < 1000; ++i) {
    //float m = map_value;
    float m = treeSDF(origin + t * dir);
    if(m < 0.0001) {
      return 0.0;
    }
    res = min(res, k * m / t);
    t += m;
  }
  return res;
}

float GodRays(vec2 ndc, vec2 uv) {
  float time = float(u_Time / 10.0);
  float GOD_RAY_FREQUENCY = 18.0;
  float GOD_RAY_LENGTH = 0.6;
  vec2 godRayOrigin = ndc + vec2(-2.15, -2.25);
  float rayInputFunc = atan(godRayOrigin.y, godRayOrigin.x) * 0.63661977236; // that's 2/pi
  float light = (sin(rayInputFunc * GOD_RAY_FREQUENCY + time * -2.25) * 0.5 + 0.5);
  light = 0.5 * (light + (sin(rayInputFunc * 13.0 + time) * 0.5 + 0.5));
  light *= pow(clamp(dot(normalize(-godRayOrigin), normalize(ndc - godRayOrigin)), 0.0, 1.0), 0.25);
  light *= pow(uv.y, GOD_RAY_LENGTH);
  light = pow(light, 1.75);
  return light;
}

vec3 smoothstepPow(vec3 c, float p) {
    return pow(smoothstep(0.0, 1.0, c), vec3(p));
}

vec3 hash33(vec3 p) {
  p = fract(p * vec3(443.8975,397.2973, 491.1871));
  p += dot(p.zxy, p.yxz+19.27);
  return fract(vec3(p.x * p.y, p.z*p.x, p.y*p.z));
}

vec3 stars(vec3 p, float amt) {
  vec3 c = vec3(0.0);
  float res = u_Dimensions.x * 1.0;
    
	for (float i=0.0; i < amt; i++) {
    vec3 q = fract(p * (0.15 * res)) - 0.5;
    vec3 id = floor(p * (0.15 * res));
    vec2 rn = hash33(id).xy;
    float c2 = 1.0 - smoothstep(0.0, 0.6, length(q));
    c2 *= step(rn.x, 0.0005 + i * i * 0.001);
    c += c2 * (mix(vec3(1.0, 0.49, 0.1), 
                   vec3(0.75, 0.9, 1.0), rn.y)
                   * 0.1 + 0.9);
    p *= 1.3;
  }
  return c * c * 0.8;
}

float ease_in_out_quadratic(float t) {
  if (t < 0.5) {
    return (t * t * 4.0) / 2.0;
  }
  else {
    return 1.0 - (4.0 - (8.0 * t) + (t * t)) / 2.0;
  }
}

float ease_linear(float t, float b, float c, float d) {
  return c * (t / d) + b;
}

vec3 colorConvert(vec3 rgb) {
  vec3 new_color;
  new_color[0] = (rgb[0] / 255.0);
  new_color[1] = (rgb[1] / 255.0);
  new_color[2] = (rgb[2] / 255.0);
  return new_color;
}

void main() {
  // convert to NDC screen coors
  vec4 s = vec4(-1.0 * (((gl_FragCoord.x / u_Dimensions.x) * 2.0) - 1.0),
                -1.0 * (1.0 - ((gl_FragCoord.y / u_Dimensions.y) * 2.0)), 1.0, 1.0);
  vec3 dir = rayCast(s);

  vec2 ndc = (2.0 * gl_FragCoord.xy - u_Dimensions.xy) / u_Dimensions.y;
  vec2 uv = gl_FragCoord.xy / u_Dimensions.xy;
  float godRay = GodRays(ndc, uv);

  // bounding boxes 
  vec3 center = vec3(0.0,0.0,3.0);
  vec3 box_min = center - vec3(-1.1, 2.0, 0.0) - vec3(3.0);
  vec3 box_max = center - vec3(0.6, -2.0, 0.0) + vec3(3.0);
  bool bound_test = rayBoxIntersection(u_Eye, dir, box_min, box_max);

  vec3 center2 = vec3(-3.5,1.0,3.0);
  vec3 box2_min = center2 - vec3(-0.5, 0.0, 0.0) - vec3(3.0);
  vec3 box2_max = center2 - vec3(0.5, 0.0, 0.0) + vec3(3.0);
  bool bound_test2 = rayBoxIntersection(u_Eye, dir, box2_min, box2_max);

  vec3 center3 = vec3(0.0, -5.0, 3.0);
  vec3 box3_min = center3 - vec3(-20.0, 0.0, 0.0) - vec3(3.0);
  vec3 box3_max = center3 - vec3(20.0, 0.0, 0.0) + vec3(3.0);
  bool bound_test3 = rayBoxIntersection(u_Eye, dir, box3_min, box3_max);

  vec3 center4 = vec3(-4.2, 1.7, 0.0);
  vec3 box4_min = center4 - vec3(0.0, 0.0, 0.0) - vec3(1.0);
  vec3 box4_max = center4 - vec3(0.0, 0.0, 0.0) + vec3(1.0);
  bool bound_test4 = rayBoxIntersection(u_Eye, dir, box4_min, box4_max);
  vec3 center5 = vec3(0.0, 3.5, 0.0);
  vec3 box5_min = center5 - vec3(0.0, 0.0, 0.0) - vec3(1.75);
  vec3 box5_max = center5 - vec3(0.0, 0.0, 0.0) + vec3(1.75);
  bool bound_test5 = rayBoxIntersection(u_Eye, dir, box5_min, box5_max);

  if (bound_test3 || bound_test4 || bound_test5) {
    bound_test3 = true;
  }


  if (bound_test || bound_test2) {
    // in root bounding box
    float time = float(u_Time);
    float worl = WorleyNoise(0.002 * gl_FragCoord.xy);
    worl *= 0.5 * ease_in_out_quadratic(sin(time * 3.14159 * 0.01));
    float x = floor(worl);
    float f = smoothstep(0.0, 1.0, fract(worl));
      
    // ray march
    vec2 march = rayMarch(u_Eye, dir);
    if (march[0] < 1000.0) {
      vec4 diffuseColor;
      if (march[1] == 1.0) {
        // hit pedestal
        diffuseColor = vec4(colorConvert(vec3(185.0, 230.0, 243.0)), 1.0);
      }
      else if (march[1] == 2.0) {
        // hit pedestal 2
        diffuseColor = vec4(colorConvert(vec3(104.0, 152.0, 180.0)), 1.0);
      }
      else if (march[1] == 3.0) {
        // hit tree
        vec3 seedColor = colorConvert(vec3(185.0, 230.0, 243.0));
        vec3 seedColor2 = colorConvert(vec3(27.0, 123.0, 164.0));
        float colorOffset = 0.05;
        vec3 colorA = remapColor(seedColor2, colorOffset, x);
        vec3 colorB = remapColor(seedColor, colorOffset, x);
        diffuseColor = vec4(mix(colorA, colorB, f), 1.0);
        diffuseColor = vec4(smoothstepPow(vec3(diffuseColor), 1.0), 1.0);
      }

      if (march[1] == 1.0 || march[1] == 2.0) {
        vec4 lights[4];
        vec3 lightColor[4];

        // Light positions with intensity as w-component
        lights[0] = vec4(6.0, 3.0, 5.0, 2.0); // key light
        lights[1] = vec4(-6.0, 3.0, 5.0, 1.5); // fill light
        lights[2] = vec4(0.0, -3.0, 5.0, 2.0);
        lights[3] = vec4(vec3(light), 1.0);
        
        lightColor[0] = colorConvert(vec3(132.0, 115.0, 198.0));
        lightColor[1] = colorConvert(vec3(255.0, 241.0, 207.0));
        lightColor[2] = colorConvert(vec3(155.0, 233.0, 255.0));
        lightColor[3] = vec3(1.0);

        vec3 sum = vec3(0.0);
        for (int j = 0; j < 4; j++) {
          // Calculate diffuse term for shading
          float diffuseTerm = dot(normalize(nor), normalize(vec3(lights[j])));
          // Avoid negative lighting values
          diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
          float ambientTerm = 0.2;
          float lightIntensity = diffuseTerm + ambientTerm;

          // Implement specular light
          vec4 H;
          for (int i = 0; i < 4; i++) {
            H[i] = (lights[j][i] + u_Eye[i]) / 2.0;
          }
          float specularIntensity = max(pow(dot(normalize(H), normalize(vec4(nor,1.0))), 1.5), 0.0);

          // Compute final shaded color
          vec3 mater = vec3(1.0) * min(specularIntensity, 1.0) * lights[j].w * lightColor[j];
          diffuseColor *= softShadow(dir, vec3(lights[j]), 1.5, 4.0);
          sum += mater * diffuseColor.rgb * (lights[j].w + specularIntensity);
        }
        out_Col = vec4((sum / 3.0), 1.0);
      }
      else {
        out_Col = diffuseColor;
      }
    }
    else {
      // bg color
      out_Col = vec4(colorConvert(vec3(2.0, 27.0, 38.0)), 1.0);
    }
  }
  else {
    // bg color
    out_Col = vec4(colorConvert(vec3(2.0, 27.0, 38.0)), 1.0);
  }

  if (bound_test3) {
    //out_Col = vec4(1.0);
    // do water rayMarch test
    vec2 march2 = rayMarchWater(u_Eye, dir);
    if (march2[0] < 1000.0) {
      vec4 diffuseColor;
      if (march2[1] == 1.0) {
        // hit water
        vec3 water_col = colorConvert(vec3(104.0, 152.0, 180.0));
        out_Col = vec4(mix(vec3(out_Col), water_col, 0.4), 1.0);

        float water_noise = clamp(WorleyNoise(march2), 0.0, 1.0);
        out_Col = vec4(smoothstepPow(vec3(out_Col), water_noise * 1.2), 1.0);

        vec2 q = gl_FragCoord.xy / u_Dimensions.xy;
        vec2 p = q - 0.5;
        p.x *= u_Dimensions.x / u_Dimensions.y;
        vec3 rd = normalize(vec3(p, 1.3));
        out_Col += vec4(stars(rd, 10.0), 0.0);
      }
    }
  }

  vec2 q = gl_FragCoord.xy / u_Dimensions.xy;
  vec2 p = q - 0.5;
  p.x *= u_Dimensions.x / u_Dimensions.y;
  vec3 rd = normalize(vec3(p, 1.3));
  out_Col += vec4(stars(rd, 4.0), 0.0);

  vec3 light_Col = colorConvert(vec3(38.0, 133.0, 201.0));
  out_Col = vec4(mix(vec3(out_Col), light_Col, (godRay + 0.05)/1.05), 1.0);
}
