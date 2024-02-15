#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

varying vec4 vertColor;
varying float Array[10];

void main() {
    gl_FragColor = vertColor;
}