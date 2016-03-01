Clazz.declarePackage ("J.renderspecial");
Clazz.load (["J.render.ShapeRenderer"], "J.renderspecial.PolyhedraRenderer", ["JU.P3", "$.V3", "JM.Atom", "JU.C"], function () {
c$ = Clazz.decorateAsClass (function () {
this.$drawEdges = 0;
this.isAll = false;
this.frontOnly = false;
this.screens3f = null;
this.scrVib = null;
this.vibs = false;
this.bsSelected = null;
this.showNumbers = false;
Clazz.instantialize (this, arguments);
}, J.renderspecial, "PolyhedraRenderer", J.render.ShapeRenderer);
Clazz.overrideMethod (c$, "render", 
function () {
var polyhedra = this.shape;
var polyhedrons = polyhedra.polyhedrons;
this.$drawEdges = polyhedra.drawEdges;
this.bsSelected = (this.vwr.getSelectionHalosEnabled () ? this.vwr.bsA () : null);
this.g3d.addRenderer (1073742182);
this.vibs = (this.ms.vibrations != null && this.tm.vibrationOn);
this.showNumbers = this.vwr.getTestFlag (3);
var needTranslucent = false;
for (var i = polyhedra.polyhedronCount; --i >= 0; ) if (polyhedrons[i].isValid && this.render1 (polyhedrons[i])) needTranslucent = true;

return needTranslucent;
});
Clazz.defineMethod (c$, "render1", 
 function (p) {
if (p.visibilityFlags == 0) return false;
var colixes = (this.shape).colixes;
var iAtom = -1;
var colix;
var scale = 1;
if (p.id == null) {
iAtom = p.centralAtom.i;
colix = (colixes == null || iAtom >= colixes.length ? 0 : colixes[iAtom]);
colix = JU.C.getColixInherited (colix, p.centralAtom.colixAtom);
} else {
colix = p.colix;
scale = p.scale;
}var needTranslucent = false;
if (JU.C.renderPass2 (colix)) {
needTranslucent = true;
} else if (!this.g3d.setC (colix)) {
return false;
}var vertices = p.vertices;
if (scale != 1) {
var v =  new Array (vertices.length);
if (scale < 0) {
var a = JU.V3.newV (p.center);
a.scale (-scale - 1);
for (var i = v.length; --i >= 0; ) {
var b = JU.V3.newV (vertices[i]);
b.add (a);
v[i] = b;
}
} else {
for (var i = v.length; --i >= 0; ) {
var a = JU.V3.newVsub (vertices[i], p.center);
a.scaleAdd2 (scale, a, p.center);
v[i] = a;
}
}vertices = v;
}if (this.screens3f == null || this.screens3f.length < vertices.length) {
this.screens3f =  new Array (vertices.length);
for (var i = vertices.length; --i >= 0; ) this.screens3f[i] =  new JU.P3 ();

}var sc = this.screens3f;
var planes = p.triangles;
for (var i = vertices.length; --i >= 0; ) {
var atom = (Clazz.instanceOf (vertices[i], JM.Atom) ? vertices[i] : null);
if (atom == null) {
this.tm.transformPtScrT3 (vertices[i], sc[i]);
} else if (atom.isVisible (this.myVisibilityFlag)) {
sc[i].set (atom.sX, atom.sY, atom.sZ);
} else if (this.vibs && atom.hasVibration ()) {
this.scrVib = this.tm.transformPtVib (atom, this.ms.vibrations[atom.i]);
sc[i].set (this.scrVib.x, this.scrVib.y, this.scrVib.z);
} else {
this.tm.transformPt3f (atom, sc[i]);
}if (this.showNumbers) {
this.g3d.setC (4);
this.g3d.drawStringNoSlab ("" + i, null, Clazz.floatToInt (sc[i].x), Clazz.floatToInt (sc[i].y), Clazz.floatToInt (sc[i].z) - 30, 0);
this.g3d.setC (colix);
}}
this.isAll = (this.$drawEdges == 1 || this.bsSelected != null);
this.frontOnly = (this.$drawEdges == 2);
var normixes = p.getNormixes ();
if (!needTranslucent || this.g3d.setC (colix)) for (var i = planes.length; --i >= 0; ) {
var pl = planes[i];
try {
if (!this.showNumbers || this.g3d.setC ((Math.round (Math.random () * 10) + 5))) this.g3d.fillTriangleTwoSided (normixes[i], sc[pl[0]], sc[pl[1]], sc[pl[2]]);
} catch (e) {
if (Clazz.exceptionOf (e, Exception)) {
System.out.println ("PolyhedraRendererError");
} else {
throw e;
}
}
}
if (this.bsSelected != null && this.bsSelected.get (iAtom)) colix = 23;
 else if (p.colixEdge != 0) colix = p.colixEdge;
if (this.g3d.setC (JU.C.getColixTranslucent3 (colix, false, 0))) for (var i = planes.length; --i >= 0; ) {
var pl = planes[i];
this.drawEdges (normixes[i], sc[pl[0]], sc[pl[1]], sc[pl[2]], -pl[3]);
}
return needTranslucent;
}, "J.shapespecial.Polyhedron");
Clazz.defineMethod (c$, "drawEdges", 
 function (normix, a, b, c, edgeMask) {
if (this.isAll || this.frontOnly && this.vwr.gdata.isDirectedTowardsCamera (normix)) {
var d = (this.g3d.isAntialiased () ? 6 : 3);
if ((edgeMask & 1) == 1) this.g3d.fillCylinderBits (3, d, a, b);
if ((edgeMask & 2) == 2) this.g3d.fillCylinderBits (3, d, b, c);
if ((edgeMask & 4) == 4) this.g3d.fillCylinderBits (3, d, a, c);
}}, "~N,JU.P3,JU.P3,JU.P3,~N");
});
