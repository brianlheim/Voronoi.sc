// Brian Heim, June 2016
// Voronoi diagram generator
// Free to use, copy, and modify.
// Written for SuperCollider 3.7.1 (Mac OS).

Voronoi : Window {
	var <points;
	var <pointBounds;
	var <colors;
	var <scaledPoints;
	var colorPoints, image;

	classvar <>printProgress = false;

	var <>drawPoints = true;
	var <>drawRegions = true;
	var <>drawPointNumbers = false;
	var <>pointRadius = 2;

	var <distFunc;

	*new {
		arg name, bounds, points, colors = nil, pointBounds = nil, distFunc = nil;
		^super.new(name, bounds, false, true).init(points, colors, pointBounds, distFunc);
	}

	init {
		arg inPoints, inColors, inPointBounds, inDistFunc;
		if(inPoints.rank==1) {inPoints = inPoints.collect(_.bubble)};
		points = inPoints;
		pointBounds = inPointBounds?this.calcPointBounds(points);
		colors = inColors?this.initColors(points);
		scaledPoints = this.scalePoints(points, pointBounds);
		distFunc = inDistFunc;
		colorPoints = Array2D.new(this.bounds.width.asInteger, this.bounds.height.asInteger);

		this.view.background_(Color.white);

		this.initDrawFunc(scaledPoints, colors);

		^this;
	}

	calcPointBounds {
		arg points;
		var xmax, xmin, ymax, ymin;
		var rmax, xcenter, ycenter;
		var pointBounds;
		var paddingFactor = 1.1;

		points.do {
			|arr|
			arr.do {
				|point|
				xmax = xmax.isNil.if {point.x} {max(xmax, point.x)};
				xmin = xmin.isNil.if {point.x} {min(xmin, point.x)};
				ymax = ymax.isNil.if {point.y} {max(ymax, point.y)};
				ymin = ymin.isNil.if {point.y} {min(ymin, point.y)};
			}
		};

		rmax = max(xmax-xmin/this.bounds.width, ymax-ymin/this.bounds.height);
		xcenter = (xmax+xmin)/2;
		ycenter = (ymax+ymin)/2;

		pointBounds = if(rmax==0) {
			Rect(xcenter-1,ycenter-1,2,2) // arbitrary, we only have one point
		} {
			var rpadded = rmax * paddingFactor;
			Rect(
				xcenter - (rpadded * this.bounds.width / 2),
				ycenter - (rpadded * this.bounds.height / 2),
				rpadded * this.bounds.width,
				rpadded * this.bounds.height
			)
		};

		//format("pointBounds: %", pointBounds).postln;
		^pointBounds;
	}

	scalePoints {
		arg points, pointBounds;

		var scaled = List[];
		points.do {
			|arr|
			scaled = scaled.add([]);
			arr.do {
				|point|
				scaled[scaled.size-1] = scaled.last add: Point(
					point.x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
					point.y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
				);
			}
		}

		^scaled;
	}

	initColors {
		arg points;

		var colors = [];

		points.do {colors = colors.add(Color.rand(0.2,0.9))};

		^colors;
	}

	distance {
		arg p1, x, y;
		distFunc.isNil.if {
			^sqrt((p1.x-x).squared + (p1.y-y).squared)
		} {
			^distFunc.value(p1.x, p1.y, x, y);
		};
	}

	colors_ {
		arg inColors;
		colors = inColors;
		this.initDrawFunc(scaledPoints, colors, false);
		^this;
	}

	distFunc_ {
		arg inDistFunc;
		distFunc = inDistFunc;
		this.initDrawFunc(scaledPoints, colors, true);
		^this;
	}

	initDrawFunc {
		arg points, colors, recalculate = true;

		// create color points
		if(recalculate) {colorPoints = this.calcColorPoints(points, colors)};


		this.drawFunc_({
			var prevcolor = -1;
			if(drawRegions) {
				this.bounds.width.asInteger.do {
					|x|
					(((x%10)==0)&&printProgress).if {x.postln};
					this.bounds.height.asInteger.do {
						|y|
						var colorindex = colorPoints.at(x,y);
						if(colorindex != prevcolor) {
							Pen.strokeColor_(colors[colorindex]);
							prevcolor = colorindex;
						};
						Pen.strokeRect(Rect(x,y,1,1));
					}
				};
			};

			if(drawPoints) {
				Pen.fillColor_(Color.black);
				points.flat.do {
					|point|
					Pen.fillOval(Rect.aboutPoint(point,pointRadius,pointRadius));
				}
			};

			if(drawPointNumbers) {
				var font = Font.default.size_(10);
				var subpoint = Point(5@5);
				if(points.maxSizeAtDepth(1) == 1) {
					points.flat.do {
						|point,i|
						Pen.stringAtPoint(i.asString,point,font);
					}
				} {

					points.do {
						|arr, i|
						arr.do {
							|point, j|
							var string = format("%,%",i,j);
							Pen.stringAtPoint(string,point,font);
						}
					}
				}
			};
		});

	}

	calcColorPoints {
		arg points;
		for(0,this.bounds.width.asInteger-1) {
			|x|
			(((x%10)==0) && printProgress).if {x.postln};
			for(0,this.bounds.height.asInteger-1) {
				|y|
				var theIndex, minDist;
				points.do {
					arg array, i;
					array.do {
						|point|
						theIndex.isNil.if {
							theIndex = i;
							minDist = this.distance(point, x, y);
						} {
							var dist = this.distance(point, x, y);
							if(dist < minDist) {
								theIndex = i;
								minDist = dist;
							}
						}
					}
				};
				colorPoints.put(x, y, theIndex);
			}
		};
		^colorPoints;
	}
}