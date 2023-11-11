PShader shader ;
PShape box ;

void setup(){
  size(600, 600, P3D);
  shader = loadShader("frag.glsl","vertex.glsl");
  
  
}

void draw(){
  background(0);
  
  color(0);
  shader(shader);
  
  translate(width/2,height/2);
  
  //shape(box);
  box(10,20,15);
  translate(30,0);
  sphere(15);
}
