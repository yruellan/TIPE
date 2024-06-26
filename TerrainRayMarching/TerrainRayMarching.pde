PShader shader;
String file = "terrain.glsl";
String file_ex = "example.glsl";

float x,y ;
float angleX, angleY ;

int scene = 0 ;
boolean bases_vectors = false ; 
boolean show_example = false ;
int variable = 0 ;

int picture_n ;
boolean take_picture ;
boolean reload_shader ;

void setup() {
  size(800, 800, P2D);
  //size(600, 600, P2D);
  noStroke();
  println("run Test.pde");
  
  shader = loadShader(file);
  
  frameRate(30);
  

  
  x = 0.0 ;
  y = -5.0 ;
  angleX = 0 ;
  angleY = 0 ;
  picture_n = 0 ;
  take_picture = false ;
  reload_shader = true ;
  
  setup_table();
  
  noLoop();
  redraw();
}

void picture(){
  float t0 = millis();
  save("../pictures/picture"+picture_n+".png");
  picture_n++ ;
  save_table();
  float t1 = millis();
  println("picture", picture_n, "(",(t1-t0)/1000.0,"s)");
  take_picture = false ;
}

void set_shader(){
  print("Reload");
  float t0 = millis();
  flush();
  shader.set("u_mouse", width/2.0+angleX, height/2.0+angleY);
  shader.set("variable", variable);

  shader.set("position", x,y);
  
  shader.set("scene",scene);
  shader.set("bases_vectors",bases_vectors);
  
  shader.set("u_resolution", float(width), float(height));
  // shader.set("sagittarius_A",loadImage("content/sagittarius_A.jpg")); 
  
  shader(shader);
  reload_shader = false ;
  float t1 = millis();
  println("\tShader Reloaded (",(t1-t0)/1000.0,"s)");
}

void draw() {
 
  if (take_picture){
    picture();
  }
  if (reload_shader){
    set_shader();
  }
  
  float t0 = millis();
  rect(0,0,width,height);
  float t1 = millis();
  println("Scene",scene,"Drawed (",(t1-t0)/1000.0,"s)");
  
}

void keyPressed(){
  float ds = 1 ;
  //println("key",keyCode,key);

  if ( keyCode >= 37 && keyCode <= 40){
    PVector spe = new PVector(0,0);
    if      (keyCode == 37){ spe = PVector.fromAngle(PI);}
    else if (keyCode == 38){ spe = PVector.fromAngle(HALF_PI);}
    else if (keyCode == 39){ spe = PVector.fromAngle(0);}
    else if (keyCode == 40){ spe = PVector.fromAngle(-HALF_PI);}
    //float t = (mouseX+angleX)/width - .5 ;
    float t = angleX/width ;
    //println("rot",t,2*t,TWO_PI*t);
    spe.rotate(-TWO_PI*t);
    x += spe.x * ds ; y += spe.y * ds ;
    shader.set("position", x,y);
  }
  else if (key == 'z' || key == 'q' ||key == 's' ||key == 'd'){
    float dtheta = 20 ;
    if      (key == 'z'){ angleY -= dtheta;}
    else if (key == 'q'){ angleX -= dtheta;}
    else if (key == 's'){ angleY += dtheta;}
    else if (key == 'd'){ angleX += dtheta;}
    
    shader.set("u_mouse", width/2.0+angleX, height/2.0+angleY);
  }
  
  else if (keyCode == 64){ // @
    bases_vectors = !bases_vectors ;
    shader.set("bases_vectors",bases_vectors);
  }
  
  else if (keyCode == 46){  // :
    variable -= 1 ;
    shader.set("variable", variable);
  }
  else if (keyCode == 47){  // =
    variable += 1 ;
    shader.set("variable", variable);
  }
  
  else if (keyCode == 92){ // `£
    show_example = !show_example ;
    println("show example :",show_example);
    if (show_example) shader = loadShader(file_ex);
    else shader = loadShader(file);
    reload_shader = true ;
  }
  else if (keyCode == 49){ // 1&
    scene = 1 ;
    shader.set("scene",scene);
  }
  else if (keyCode == 50){ // 2é
    scene = 2 ;
    shader.set("scene",scene);
  }
  else if (keyCode == 51){ // 3"
    scene = 3 ;
    shader.set("scene",scene);
  }
  else if (keyCode == 52){ // 4'
    scene = 4 ;
    shader.set("scene",scene);
  }
  else if (key == 'l'){
    scene -= 1 ;
    shader.set("scene",scene);
  }
  else if (key == 'm'){
    scene += 1 ;
    shader.set("scene",scene);
  }

  else if (keyCode == 10){ // ENTER
    take_picture = true ;
  }
  else if (keyCode == 27){ // ESC
    exit();
  }
  else if (keyCode == 84){ // t
    println("Time : ",millis()/1000.0);
  }
  else if (keyCode == 32){ // SPACE
    shader = loadShader(file);
    reload_shader = true ;
  } else {
    println("Key : ",key,": ",keyCode);
  }
  redraw();
}
