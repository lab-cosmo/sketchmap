Clazz.declarePackage ("JS");
Clazz.load (null, "JS.SmilesStereo", ["java.lang.Float", "java.util.Arrays", "JU.AU", "$.Measure", "$.PT", "$.T3", "$.V3", "JS.InvalidSmilesException", "$.PolyhedronStereoSorter", "$.SmilesAromatic", "$.SmilesAtom", "$.SmilesParser", "JU.Escape", "$.Logger"], function () {
c$ = Clazz.decorateAsClass (function () {
this.chiralClass = -2147483648;
this.chiralOrder = -2147483648;
this.atomCount = 0;
this.details = null;
this.search = null;
this.jmolAtoms = null;
this.directives = null;
this.polyhedralOrders = null;
this.isNot = false;
this.sorter = null;
Clazz.instantialize (this, arguments);
}, JS, "SmilesStereo");
c$.getChiralityClass = Clazz.defineMethod (c$, "getChiralityClass", 
 function (xx) {
return Clazz.doubleToInt (("0;PH;AL;33;TH;TP;OH;77;SP;".indexOf (xx) + 1) / 3);
}, "~S");
c$.newStereo = Clazz.defineMethod (c$, "newStereo", 
function (stereo) {
return (stereo == null ?  new JS.SmilesStereo (0, 0, 0, null, null) :  new JS.SmilesStereo (stereo.chiralClass, stereo.chiralOrder, stereo.atomCount, stereo.details, stereo.directives));
}, "JS.SmilesStereo");
Clazz.makeConstructor (c$, 
function (chiralClass, chiralOrder, atomCount, details, directives) {
this.chiralClass = chiralClass;
this.chiralOrder = chiralOrder;
this.atomCount = atomCount;
this.details = details;
this.directives = directives;
if (chiralClass == 1) this.getPolyhedralOrders ();
}, "~N,~N,~N,~S,~S");
Clazz.defineMethod (c$, "getPolyhedralOrders", 
 function () {
var po = this.polyhedralOrders = JU.AU.newInt2 (this.atomCount);
if (this.details == null) return;
var temp =  Clazz.newIntArray (this.details.length, 0);
var ret =  Clazz.newIntArray (1, 0);
var msg = null;
var pt = 0;
var s = this.details + "/";
var n = 0;
var len = s.length;
var index = 0;
var atomPt = 0;
do {
var ch = s.charAt (index);
switch (ch) {
case '!':
this.isNot = true;
index++;
break;
case '/':
case '.':
if ((pt = atomPt) >= this.atomCount) {
msg = "Too many descriptors";
break;
}var a = po[atomPt] =  Clazz.newIntArray (n, 0);
for (; --n >= 0; ) a[n] = temp[n];

n = 0;
if (JU.Logger.debugging) JU.Logger.info (JU.PT.toJSON ("@PH" + this.atomCount + "[" + atomPt + "]", a));
if (ch == '/') index = 2147483647;
 else index++;
atomPt++;
break;
default:
index = JS.SmilesParser.getRingNumber (s, index, ch, ret);
pt = temp[n++] = ret[0] - 1;
if (pt == atomPt) msg = "Atom cannot connect to itself";
 else if (pt < 0 || pt >= this.atomCount) msg = "Connection number outside of range (1-" + this.atomCount + ")";
 else if (n >= this.atomCount) msg = "Too many connections indicated";
}
if (msg != null) {
msg += ": " + s.substring (0, index) + "<<";
throw  new JS.InvalidSmilesException (msg);
}} while (index < len);
});
Clazz.defineMethod (c$, "getChiralClass", 
function () {
return this.chiralClass;
});
Clazz.defineMethod (c$, "setChiralClass", 
function (chiralClass) {
this.chiralClass = chiralClass;
}, "~N");
Clazz.defineMethod (c$, "getChiralOrder", 
function () {
return this.chiralOrder;
});
Clazz.defineMethod (c$, "setChiralOrder", 
function (chiralOrder) {
this.chiralOrder = chiralOrder;
}, "~N");
Clazz.defineMethod (c$, "fixStereo", 
function (sAtom) {
var nBonds = Math.max (sAtom.missingHydrogenCount, 0) + sAtom.getBondCount ();
switch (this.chiralClass) {
case 0:
switch (nBonds) {
case 2:
this.chiralClass = 2;
break;
case 3:
this.chiralClass = 3;
break;
case 4:
case 5:
case 6:
this.chiralClass = nBonds;
break;
}
break;
case 8:
if (nBonds != 4) sAtom.stereo = null;
break;
case 1:
if (nBonds != 0 && nBonds != this.atomCount) sAtom.stereo = null;
break;
case 2:
case 6:
case 4:
case 5:
if (nBonds != this.chiralClass) sAtom.stereo = null;
break;
}
if (sAtom.stereo == null) throw  new JS.InvalidSmilesException ("Incorrect number of bonds for stereochemistry descriptor");
}, "JS.SmilesAtom");
Clazz.defineMethod (c$, "setSmilesCoordinates", 
 function (atom, sAtom, sAtom2, cAtoms) {
if (atom.stereo == null) return false;
var chiralClass = atom.stereo.chiralClass;
var chiralOrder = atom.stereo.chiralOrder;
var a2 = (chiralClass == 2 || chiralClass == 3 ? a2 = this.jmolAtoms[sAtom2.getMatchingAtomIndex ()] : null);
atom.set (0, 0, 0);
atom = this.jmolAtoms[sAtom.getMatchingAtomIndex ()];
atom.set (0, 0, 0);
var map = this.search.getMappedAtoms (atom, a2, cAtoms);
switch (chiralClass) {
case 1:
break;
case 2:
case 4:
if (chiralOrder == 2) {
var i = map[0];
map[0] = map[1];
map[1] = i;
}cAtoms[map[0]].set (0, 0, 1);
cAtoms[map[1]].set (1, 0, -1);
cAtoms[map[2]].set (0, 1, -1);
cAtoms[map[3]].set (-1, -1, -1);
break;
case 8:
switch (chiralOrder) {
case 1:
cAtoms[map[0]].set (1, 0, 0);
cAtoms[map[1]].set (0, 1, 0);
cAtoms[map[2]].set (-1, 0, 0);
cAtoms[map[3]].set (0, -1, 0);
break;
case 2:
cAtoms[map[0]].set (1, 0, 0);
cAtoms[map[1]].set (-1, 0, 0);
cAtoms[map[2]].set (0, 1, 0);
cAtoms[map[3]].set (0, -1, 0);
break;
case 3:
cAtoms[map[0]].set (1, 0, 0);
cAtoms[map[1]].set (0, 1, 0);
cAtoms[map[2]].set (0, -1, 0);
cAtoms[map[3]].set (-1, 0, 0);
break;
}
break;
case 5:
case 6:
var n = map.length;
if (chiralOrder == 2) {
var i = map[0];
map[0] = map[n - 1];
map[n - 1] = i;
}cAtoms[map[0]].set (0, 0, 1);
cAtoms[map[n - 1]].set (0, 0, -1);
cAtoms[map[1]].set (1, 0, 0);
cAtoms[map[2]].set (0, 1, 0);
cAtoms[map[3]].set (-1, 0, 0);
if (n == 6) cAtoms[map[4]].set (0, -1, 0);
break;
}
return true;
}, "JS.SmilesAtom,JS.SmilesAtom,JS.SmilesAtom,~A");
Clazz.defineMethod (c$, "getX", 
 function (sAtom, jn, pt, haveCoordinates, needHSwitch) {
var atom = this.getJmolAtom (sAtom.getMatchingAtomIndex ());
var doSwitch = sAtom.isFirst || pt == 3;
if (haveCoordinates) {
if (this.search.isSmarts) {
var b = atom.getEdges ();
for (var i = 0; i < b.length; i++) {
if (b[i].getCovalentOrder () == 2) continue;
var a = this.jmolAtoms[atom.getBondedAtomIndex (i)];
if (a === jn[pt - 1]) continue;
jn[pt] = a;
break;
}
}if (jn[pt] == null) {
var v =  new JU.V3 ();
var n = 0;
for (var i = 0; i < 4; i++) {
if (jn[i] == null) continue;
n++;
v.sub (jn[i]);
}
if (v.length () == 0) {
v.setT ((jn[4]));
doSwitch = false;
} else {
v.scaleAdd2 (n + 1, this.getJmolAtom (sAtom.getMatchingAtomIndex ()), v);
doSwitch = this.search.isSmilesFind || doSwitch;
}jn[pt] =  new JS.SmilesAtom ().setIndex (-1);
(jn[pt]).setT (v);
}}if (jn[pt] == null) {
jn[pt] = this.search.getHydrogens (atom, null);
if (needHSwitch) doSwitch = true;
}if (jn[pt] != null && doSwitch) {
var a = jn[pt];
jn[pt] = jn[pt - 1];
jn[pt - 1] = a;
}}, "JS.SmilesAtom,~A,~N,~B,~B");
Clazz.defineMethod (c$, "getJmolAtom", 
 function (i) {
return (i < 0 || i >= this.jmolAtoms.length ? null : this.jmolAtoms[i]);
}, "~N");
Clazz.defineMethod (c$, "sortBondsByStereo", 
function (atom, atomPrev, ref, bonds, vTemp) {
if (bonds.length < 2 || !(Clazz.instanceOf (atom, JU.T3))) return;
if (atomPrev == null) atomPrev = bonds[0].getOtherAtomNode (atom);
var aTemp =  Clazz.newArray (bonds.length, 0, null);
if (this.sorter == null) this.sorter =  new JS.PolyhedronStereoSorter ();
vTemp.sub2 (atomPrev, ref);
this.sorter.setRef (vTemp);
for (var i = bonds.length; --i >= 0; ) {
var a = bonds[i].getOtherAtomNode (atom);
var f = (a === atomPrev ? 0 : this.sorter.isAligned (a, ref, atomPrev) ? -999 : JU.Measure.computeTorsion (atom, atomPrev, ref, a, true));
if (bonds.length > 2) f += 360;
aTemp[i] =  Clazz.newArray (-1, [bonds[i], Float.$valueOf (f), a]);
}
java.util.Arrays.sort (aTemp, this.sorter);
if (JU.Logger.debugging) JU.Logger.info (JU.Escape.e (aTemp));
for (var i = bonds.length; --i >= 0; ) bonds[i] = aTemp[i][0];

}, "JU.Node,JU.Node,JU.T3,~A,JU.V3");
Clazz.defineMethod (c$, "checkStereoChemistry", 
function (smilesSearch, v) {
this.search = smilesSearch;
this.jmolAtoms = this.search.jmolAtoms;
var isSmilesFind = smilesSearch.isSmilesFind;
var invertStereochemistry = smilesSearch.invertStereochemistry;
if (JU.Logger.debugging) JU.Logger.debug ("checking sstereochemistry...");
for (var i = 0; i < smilesSearch.ac; i++) {
var sAtom = smilesSearch.patternAtoms[i];
if (sAtom.stereo == null) continue;
var isNot = (sAtom.not != invertStereochemistry);
var atom1 = null;
var atom2 = null;
var atom3 = null;
var atom4 = null;
var atom5 = null;
var atom6 = null;
var sAtom1 = null;
var sAtom2 = null;
var sAtom0 = null;
var jn;
var atom0 = sAtom.getMatchingAtom ();
if (isSmilesFind) sAtom0 = atom0;
var nH = Math.max (sAtom.missingHydrogenCount, 0);
var order = sAtom.stereo.chiralOrder;
var chiralClass = sAtom.stereo.chiralClass;
if (isSmilesFind && sAtom0.getChiralClass () != chiralClass) return false;
if (JU.Logger.debugging) JU.Logger.debug ("...type " + chiralClass + " for pattern atom " + sAtom + " " + atom0);
switch (chiralClass) {
case 1:
if (sAtom.stereo.isNot) isNot = !isNot;
if (nH > 1 || sAtom.bondCount == 0) continue;
if (isSmilesFind) {
continue;
}var bonds = sAtom.bonds;
var jHpt = -1;
if (nH == 1) {
jHpt = (sAtom.isFirst ? 0 : 1);
if (sAtom.getBondCount () != 3) return false;
v.vA.set (0, 0, 0);
for (var j = 0; j < 3; j++) v.vA.add (bonds[j].getOtherAtom (sAtom0).getMatchingAtom ());

v.vA.scale (0.3333);
v.vA.sub2 (atom0, v.vA);
v.vA.add (atom0);
}var po = sAtom.stereo.polyhedralOrders;
var pt;
for (var j = po.length; --j >= 0; ) {
var orders = po[j];
if (orders == null || orders.length < 2) continue;
pt = (j > jHpt ? j - nH : j);
var ta1 = (j == jHpt ? v.vA : bonds[pt].getOtherAtom (sAtom).getMatchingAtom ());
var flast = (isNot ? 3.4028235E38 : 0);
var ta2 = null;
for (var k = 0; k < orders.length; k++) {
pt = orders[k];
var ta3;
if (pt == jHpt) {
ta3 = v.vA;
} else {
if (pt > jHpt) pt--;
ta3 = bonds[pt].getOtherAtom (sAtom).getMatchingAtom ();
}if (k == 0) {
ta2 = ta3;
continue;
}var f = JU.Measure.computeTorsion (ta3, ta1, atom0, ta2, true);
if (Float.isNaN (f)) f = 180;
if (orders.length == 2) return ((f < 0) != isNot);
if (f < 0) f += 360;
if ((f < flast) != isNot) return false;
flast = f;
}
}
continue;
case 2:
var isAllene = true;
if (isAllene) {
sAtom1 = sAtom.getBond (0).getOtherAtom (sAtom);
sAtom2 = sAtom.getBond (1).getOtherAtom (sAtom);
if (sAtom1 == null || sAtom2 == null) continue;
var sAtom1a = sAtom;
var sAtom2a = sAtom;
while (sAtom1.getBondCount () == 2 && sAtom2.getBondCount () == 2 && sAtom1.getValence () == 4 && sAtom2.getValence () == 4) {
var b = sAtom1.getBondNotTo (sAtom1a, true);
sAtom1a = sAtom1;
sAtom1 = b.getOtherAtom (sAtom1);
b = sAtom2.getBondNotTo (sAtom2a, true);
sAtom2a = sAtom2;
sAtom2 = b.getOtherAtom (sAtom2);
}
sAtom = sAtom1;
}jn =  new Array (6);
jn[4] =  new JS.SmilesAtom ().setIndex (604);
var nBonds = sAtom.getBondCount ();
for (var k = 0; k < nBonds; k++) {
sAtom1 = sAtom.bonds[k].getOtherAtom (sAtom);
if (sAtom.bonds[k].matchingBond.getCovalentOrder () == 2) {
if (sAtom2 == null) sAtom2 = sAtom1;
} else if (jn[0] == null) {
jn[0] = sAtom1.getMatchingAtom ();
} else {
jn[1] = sAtom1.getMatchingAtom ();
}}
if (sAtom2 == null) continue;
nBonds = sAtom2.getBondCount ();
if (nBonds < 2 || nBonds > 3) continue;
for (var k = 0; k < nBonds; k++) {
sAtom1 = sAtom2.bonds[k].getOtherAtom (sAtom2);
if (sAtom2.bonds[k].matchingBond.getCovalentOrder () == 2) {
} else if (jn[2] == null) {
jn[2] = sAtom1.getMatchingAtom ();
} else {
jn[3] = sAtom1.getMatchingAtom ();
}}
if (isSmilesFind) {
if (jn[1] == null) this.getX (sAtom, jn, 1, false, isAllene);
if (jn[3] == null) this.getX (sAtom2, jn, 3, false, false);
if (!this.setSmilesCoordinates (sAtom0, sAtom, sAtom2, jn)) return false;
}if (jn[1] == null) this.getX (sAtom, jn, 1, true, false);
if (jn[3] == null) this.getX (sAtom2, jn, 3, true, false);
if (!JS.SmilesStereo.checkStereochemistryAll (sAtom.not != invertStereochemistry, atom0, chiralClass, order, jn[0], jn[1], jn[2], jn[3], null, null, v)) return false;
continue;
case 3:
case 4:
case 8:
case 5:
case 6:
atom1 = this.getJmolAtom (sAtom.getMatchingBondedAtom (0));
switch (nH) {
case 0:
atom2 = this.getJmolAtom (sAtom.getMatchingBondedAtom (1));
break;
case 1:
atom2 = smilesSearch.getHydrogens (sAtom.getMatchingAtom (), null);
if (sAtom.isFirst) {
var a = atom2;
atom2 = atom1;
atom1 = a;
}break;
default:
continue;
}
atom3 = this.getJmolAtom (sAtom.getMatchingBondedAtom (2 - nH));
atom4 = this.getJmolAtom (sAtom.getMatchingBondedAtom (3 - nH));
atom5 = this.getJmolAtom (sAtom.getMatchingBondedAtom (4 - nH));
atom6 = this.getJmolAtom (sAtom.getMatchingBondedAtom (5 - nH));
if (isSmilesFind && !this.setSmilesCoordinates (sAtom0, sAtom, sAtom2,  Clazz.newArray (-1, [atom1, atom2, atom3, atom4, atom5, atom6]))) return false;
if (!JS.SmilesStereo.checkStereochemistryAll (isNot, atom0, chiralClass, order, atom1, atom2, atom3, atom4, atom5, atom6, v)) return false;
continue;
}
}
return true;
}, "JS.SmilesSearch,JS.VTemp");
c$.getStereoFlag = Clazz.defineMethod (c$, "getStereoFlag", 
function (atom0, atoms, nAtoms, v) {
var atom1 = atoms[0];
var atom2 = atoms[1];
var atom3 = atoms[2];
var atom4 = atoms[3];
var atom5 = atoms[4];
var atom6 = atoms[5];
var chiralClass = 4;
switch (nAtoms) {
default:
case 5:
case 6:
return (JS.SmilesStereo.checkStereochemistryAll (false, atom0, chiralClass, 1, atom1, atom2, atom3, atom4, atom5, atom6, v) ? "@" : "@@");
case 2:
case 4:
if (atom3 == null || atom4 == null) return "";
var d = JS.SmilesAromatic.getNormalThroughPoints (atom1, atom2, atom3, v.vTemp, v.vA, v.vB);
if (Math.abs (JS.SmilesStereo.distanceToPlane (v.vTemp, d, atom4)) < 0.2) {
chiralClass = 8;
if (JS.SmilesStereo.checkStereochemistryAll (false, atom0, chiralClass, 1, atom1, atom2, atom3, atom4, atom5, atom6, v)) return "@SP1";
if (JS.SmilesStereo.checkStereochemistryAll (false, atom0, chiralClass, 2, atom1, atom2, atom3, atom4, atom5, atom6, v)) return "@SP2";
if (JS.SmilesStereo.checkStereochemistryAll (false, atom0, chiralClass, 3, atom1, atom2, atom3, atom4, atom5, atom6, v)) return "@SP3";
} else {
return (JS.SmilesStereo.checkStereochemistryAll (false, atom0, chiralClass, 1, atom1, atom2, atom3, atom4, atom5, atom6, v) ? "@" : "@@");
}}
return "";
}, "JU.Node,~A,~N,JS.VTemp");
c$.checkStereochemistryAll = Clazz.defineMethod (c$, "checkStereochemistryAll", 
function (isNot, atom0, chiralClass, order, atom1, atom2, atom3, atom4, atom5, atom6, v) {
switch (chiralClass) {
default:
return true;
case 1:
return true;
case 3:
return (isNot == (JS.SmilesStereo.getHandedness (atom2, atom3, atom0, atom1, v) != order));
case 2:
case 4:
return (isNot == (JS.SmilesStereo.getHandedness (atom2, atom3, atom4, atom1, v) != order));
case 5:
return (isNot == (!JS.SmilesStereo.isDiaxial (atom0, atom0, atom5, atom1, v, -0.95) || JS.SmilesStereo.getHandedness (atom2, atom3, atom4, atom1, v) != order));
case 6:
if (isNot != (!JS.SmilesStereo.isDiaxial (atom0, atom0, atom6, atom1, v, -0.95))) return false;
JS.SmilesStereo.getPlaneNormals (atom2, atom3, atom4, atom5, v);
if (isNot != (v.vNorm1.dot (v.vNorm2) < 0 || v.vNorm2.dot (v.vNorm3) < 0)) return false;
v.vNorm2.sub2 (atom0, atom1);
return (isNot == ((v.vNorm1.dot (v.vNorm2) < 0 ? 2 : 1) == order));
case 8:
JS.SmilesStereo.getPlaneNormals (atom1, atom2, atom3, atom4, v);
return (v.vNorm1.dot (v.vNorm2) < 0 ? isNot == (order != 3) : v.vNorm2.dot (v.vNorm3) < 0 ? isNot == (order != 2) : isNot == (order != 1));
}
}, "~B,JU.Node,~N,~N,JU.Node,JU.Node,JU.Node,JU.Node,JU.Node,JU.Node,JS.VTemp");
c$.isDiaxial = Clazz.defineMethod (c$, "isDiaxial", 
function (atomA, atomB, atom1, atom2, v, f) {
v.vA.sub2 (atomA, atom1);
v.vB.sub2 (atomB, atom2);
v.vA.normalize ();
v.vB.normalize ();
return (v.vA.dot (v.vB) < f);
}, "JU.Node,JU.Node,JU.Node,JU.Node,JS.VTemp,~N");
c$.getHandedness = Clazz.defineMethod (c$, "getHandedness", 
 function (a, b, c, pt, v) {
var d = JS.SmilesAromatic.getNormalThroughPoints (a, b, c, v.vTemp, v.vA, v.vB);
return (JS.SmilesStereo.distanceToPlane (v.vTemp, d, pt) > 0 ? 1 : 2);
}, "JU.Node,JU.Node,JU.Node,JU.Node,JS.VTemp");
c$.getPlaneNormals = Clazz.defineMethod (c$, "getPlaneNormals", 
 function (atom1, atom2, atom3, atom4, v) {
JS.SmilesAromatic.getNormalThroughPoints (atom1, atom2, atom3, v.vNorm1, v.vTemp1, v.vTemp2);
JS.SmilesAromatic.getNormalThroughPoints (atom2, atom3, atom4, v.vNorm2, v.vTemp1, v.vTemp2);
JS.SmilesAromatic.getNormalThroughPoints (atom3, atom4, atom1, v.vNorm3, v.vTemp1, v.vTemp2);
}, "JU.Node,JU.Node,JU.Node,JU.Node,JS.VTemp");
c$.distanceToPlane = Clazz.defineMethod (c$, "distanceToPlane", 
function (norm, w, pt) {
return (norm == null ? NaN : (norm.x * pt.x + norm.y * pt.y + norm.z * pt.z + w) / Math.sqrt (norm.x * norm.x + norm.y * norm.y + norm.z * norm.z));
}, "JU.V3,~N,JU.P3");
c$.checkChirality = Clazz.defineMethod (c$, "checkChirality", 
function (pattern, index, newAtom) {
var stereoClass = 0;
var order = -2147483648;
var len = pattern.length;
var details = null;
var directives = null;
var atomCount = 0;
var ch;
stereoClass = 0;
order = 1;
var isPoly = false;
if (++index < len) {
switch (ch = pattern.charAt (index)) {
case '@':
order = 2;
index++;
break;
case 'H':
break;
case 'P':
isPoly = true;
case 'A':
case 'O':
case 'S':
case 'T':
stereoClass = (index + 1 < len ? JS.SmilesStereo.getChiralityClass (pattern.substring (index, index + 2)) : -1);
index += 2;
break;
default:
order = (JU.PT.isDigit (ch) ? 1 : -1);
}
var pt = index;
if (order == 1 || isPoly) {
while (pt < len && JU.PT.isDigit (pattern.charAt (pt))) pt++;

if (pt > index) {
try {
var n = Integer.parseInt (pattern.substring (index, pt));
if (isPoly) {
atomCount = n;
if (pt < len && pattern.charAt (pt) == '(') {
details = JS.SmilesParser.getSubPattern (pattern, pt, '(');
pt += details.length + 2;
}if (pt < len && pattern.charAt (pt) == '/') {
directives = JS.SmilesParser.getSubPattern (pattern, pt, '/');
pt += directives.length + 2;
}} else {
order = n;
}} catch (e) {
if (Clazz.exceptionOf (e, NumberFormatException)) {
order = -1;
} else {
throw e;
}
}
index = pt;
}}if (order < 1 || stereoClass < 0) throw  new JS.InvalidSmilesException ("Invalid stereochemistry descriptor");
}newAtom.stereo =  new JS.SmilesStereo (stereoClass, order, atomCount, details, directives);
if (JS.SmilesParser.getChar (pattern, index) == '?') {
JU.Logger.info ("Ignoring '?' in stereochemistry");
index++;
}return index;
}, "~S,~N,JS.SmilesAtom");
Clazz.defineStatics (c$,
"STEREOCHEMISTRY_SQUARE_PLANAR", 8,
"STEREOCHEMISTRY_OCTAHEDRAL", 6,
"STEREOCHEMISTRY_TRIGONAL_BIPYRAMIDAL", 5,
"STEREOCHEMISTRY_TETRAHEDRAL", 4,
"STEREOCHEMISTRY_TRIGONAL_PYRAMIDAL", 3,
"STEREOCHEMISTRY_ALLENE", 2,
"STEREOCHEMISTRY_POLYHEDRAL", 1,
"STEREOCHEMISTRY_DEFAULT", 0);
});
