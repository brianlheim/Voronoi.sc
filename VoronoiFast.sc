VoronoiFast : Window {
	var <points, <perimeters;
	var allPerimeters;
	var <pointBounds;
	var boundaryLines;
	var <>errorAllowance;


	*new {
		arg name, bounds, points, pointBounds = nil;
		^super.new(name, bounds, false, false).init(points, pointBounds);
	}

	log {
		arg prevTime, string;
		var newTime = Process.elapsedTime;
		"%\n\tcurrent time: %\tdelta time: %".format(string, newTime, (newTime - prevTime).round(0.01)).postln;
		^newTime;
	}

	init {
		arg inPoints, inPointBounds;

		var time = this.log(Process.elapsedTime, "starting timer");
		// error checking on points
		if(inPoints.isEmpty) {Error("Points list is empty").throw};
		// check for duplicate points
		if(inPoints.size >= 2) {
			for(0,inPoints.size-2, {
				|i|
				for(i+1,inPoints.size-1, {
					|j|
					if(inPoints[i] == inPoints[j]) {
						Error("Duplicate point in points list").throw
					};
				});
			});
		};
		// time = this.log(time, "error checking");

		points = inPoints.collect({|p| p.x.asFloat@p.y.asFloat});
		perimeters = List[];
		pointBounds = inPointBounds?this.calcPointBounds(points);
		// time = this.log(time, "calc point bounds");
		this.errorAllowance = this.calcErrorAllowance(pointBounds);
		// time = this.log(time, "calc error allowance");
		boundaryLines = this.initBoundaryLines(pointBounds);
		time = this.log(time, "init boundary lines");

		time = this.log(time, "points loop starting");
		points.do {
			|point, i|
			var lines, intersections, perimeter;
			lines = this.calcLines(point, points);
			// time = this.log(time, "lines %".format(i));
			intersections = this.calcIntersections(point, lines);
			// time = this.log(time, "intersect %".format(i));
			perimeter = this.calcPerimeter(point, lines, intersections);
			// time = this.log(time, "perimeter %".format(i));
			perimeters = perimeters.add(perimeter);
		};
		time = this.log(time, "finished points loop");

		allPerimeters = this.calcUniquePerimeters(perimeters);
		time = this.log(time, "calc unique perimeters");

		this.view.background_(Color.white);

		this.initDrawFunc(points, perimeters, allPerimeters);
		time = this.log(time, "init draw func");

		this.refresh();
		time = this.log(time, "refresh");
		this.front();
		^this;
	}

	calcPointBounds {
		arg points;
		var xmax, xmin, ymax, ymin;
		var rmax, xcenter, ycenter;
		var pointBounds;
		var paddingFactor = 1.05;

		points.do {
			|point|
			xmax = xmax.isNil.if {point.x} {max(xmax, point.x)};
			xmin = xmin.isNil.if {point.x} {min(xmin, point.x)};
			ymax = ymax.isNil.if {point.y} {max(ymax, point.y)};
			ymin = ymin.isNil.if {point.y} {min(ymin, point.y)};
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

		format("pointBounds: %", pointBounds).postln;
		^pointBounds;
	}

	calcErrorAllowance {
		arg pointBounds;

		var errorAllowance = max(pointBounds.width, pointBounds.height) * 1e-10;
		format("errorAllowance: %", errorAllowance).postln;
		^errorAllowance;
	}

	initBoundaryLines {
		arg pointBounds;
		var left = [1,0,pointBounds.left.neg];
		var right = [1,0,pointBounds.right.neg];
		var top = [0,1,pointBounds.top.neg];
		var bottom = [0,1,pointBounds.bottom.neg];
		^List[left, right, top, bottom];
	}

	calcLines {
		arg point, points;
		var lines = boundaryLines.copy;

		points.do {
			|point2|
			if(point != point2) {
				// normally point.y & point2.y are switched,
				// but we want the perpendicular bisector
				var midpoint = Point(point.x+point2.x/2, point.y+point2.y/2);
				var pointA = (point-midpoint).rotate(0.5pi)+midpoint;
				var pointB = (point2-midpoint).rotate(0.5pi)+midpoint;
				var line = [
					this.thresh2(pointA.y - pointB.y),
					this.thresh2(pointB.x - pointA.x),
					this.thresh2((pointA.x*pointB.y)-(pointB.x*pointA.y))
				];
				lines = lines.add(line);
				/*format("for points (%,%), (%,%), bisector is (%,%,%)",
				point.x,point.y,point2.x,point2.y,line[0],line[1],line[2]
				).postln;*/
			}
		};

		^lines;
	}

	calcIntersections {
		arg point, lines;
		var intersections = List[];
		var externalLines = List[];

		for(0,lines.size-2) {
			|i|
			for(i+1,lines.size-1) {
				|j|
				var lineA, lineB, a1, a2, b1, b2, c1, c2, intersect;
				lineA = lines[i];
				lineB = lines[j];
				a1 = lineA[0];
				b1 = lineA[1];
				c1 = lineA[2];
				a2 = lineB[0];
				b2 = lineB[1];
				c2 = lineB[2];

				intersect = case
				{(a1==0)&&(a2==0)} {nil} // both horizontal
				{(b1==0)&&(b2==0)} {nil} // both vertical
				{(a1==0)&&(b2==0)} {Point(c2.neg/a2, c1.neg/b1)} // perp
				{(a2==0)&&(b1==0)} {Point(c1.neg/a1, c2.neg/b2)} // perp
				{b1==0} {var x = c1.neg/a1; Point(x, (a2*x+c2)/b2.neg)}
				{b2==0} {var x = c2.neg/a2; Point(x, (a1*x*c1)/b1.neg)}
				{((a1/b1)-(a2/b2))==0} {nil} // parallel
				{
					var x = ((c2/b2)-(c1/b1))/((a1/b1)-(a2/b2));
					var y = (c1+(a1*x))/b1.neg;
					Point(x,y);
				};

				/*format("for lines (%,%,%), (%,%,%), intersect is %",
				a1,b1,c1,a2,b2,c2,
				intersect.isNil.if {"none"} {
				"(%,%)".format(intersect.x, intersect.y)
				}).postln;*/

				intersect = if(intersect.notNil) {
					if(this.pointBounds.contains(intersect)) {intersect} {nil}
				} {nil};

				if(intersect.notNil) {intersections = intersections.add(intersect)};
			}
		};

		// lines.removeAll(externalLines);
		//lines.size.postln;

		// format("valid intersection points: %", intersections.size).postln;

		^intersections;
	}

	calcPerimeter {
		arg point, lines, intersections;

		var perimeter = intersections;

		// lines.do(_.postln);
		lines.do {
			|line|
			var a=line[0], b=line[1], c=line[2];
			var dir, pL, pvector0, diff;
			dir = Point(b, a.neg);
			dir = dir / dir.dist(Point(0,0));

			pL = (a==0).if {
				//a=0
				Point(0,c/b.neg)
			} {
				Point(c/a.neg,0)
			};

			diff = pL - point;

			pvector0 = point - pL + (((diff.x*dir.x)+(diff.y*dir.y))*dir);
			pvector0.x = this.thresh2(pvector0.x);
			pvector0.y = this.thresh2(pvector0.y);
			// p - pL + ((p - pL) dot dir)*dir
			//"line (%,%,%), target point (%,%) has perpendicular vector (%,%)".format(
			//	a,b,c,point.x,point.y,pvector0.x,pvector0.y).postln;

			perimeter = perimeter.reject({
				|pointB|
				var diffB, pvectorB, isWrongSide;
				diffB = pL - pointB;
				pvectorB = pointB - pL + (((diffB.x*dir.x)+(diffB.y*dir.y))*dir);
				pvectorB.x = this.thresh2(pvectorB.x);
				pvectorB.y = this.thresh2(pvectorB.y);
				//"line (%,%,%), candidate point (%,%) has perpendicular vector(%,%)".format(
				//	a,b,c,pointB.x,pointB.y,pvectorB.x,pvectorB.y).postln;
				isWrongSide = case
				{(pvectorB.x == 0) && (pvectorB.y == 0)} {false}
				{(pvectorB.x.sign != pvector0.x.sign) || (pvectorB.y.sign != pvector0.y.sign)} {true}
				{false};

				//"for vectors (%,%) and (%,%), result is: % side".format(
				//	pvector0.x,pvector0.y,pvectorB.x,pvectorB.y,isWrongSide.if {"wrong"} {"right"}).postln;

				isWrongSide;
			});
			//perimeter.size.postln;
		};

		//"perimeter for point (%,%) has perimeter of size %".format(point.x.round(0.01),point.y.round(0.01),perimeter.size).postln;

		perimeter = this.giftwrap(perimeter);
		//perimeter.collect(_.round(0.01)).postln;

		^perimeter;
	}

	calcUniquePerimeters {
		arg perimeters;
		^perimeters.reduce('++');
	}

	initDrawFunc {
		arg points, perimeters, allPerimeters;
		this.drawFunc_({
			this.drawFields(perimeters);
			this.drawPerimeters(allPerimeters);
			this.drawGrid();
			this.drawPoints(points);
		});
	}

	drawPoints {
		arg points;
		var radius = 2;
		Pen.fillColor_(Color.black);
		points.do {
			|point|
			var scaledPoint = Point(
				point.x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
				point.y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
			);
			//format("drawing point (%, %)", scaledPoint.x, scaledPoint.y).postln;
			Pen.fillOval(Rect.aboutPoint(scaledPoint,radius,radius));
		}
	}

	drawGrid {

	}

	drawFields {
		arg perimeters;
		var radius = 5;


		perimeters.do {
			|perimeter,i|
			Pen.fillColor_(Color.rand(0.6,1));
			//"drawing perimeter %".format(i).postln;
			perimeter.do {
				|pair, i|
				var scaledPoint1 = Point(
					pair[0].x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
					pair[0].y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
				);
				var scaledPoint2 = Point(
					pair[1].x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
					pair[1].y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
				);
				if(i == 0) {Pen.moveTo(scaledPoint1)};
				Pen.lineTo(scaledPoint2);
				//[scaledPoint1, scaledPoint2].postln;
			};
			Pen.fill;
		}
	}

	drawPerimeters {
		arg perimeters;
		Pen.strokeColor_(Color.black);
		perimeters.do {
			|pair|
			var scaledPoint1 = Point(
				pair[0].x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
				pair[0].y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
			);
			var scaledPoint2 = Point(
				pair[1].x.linlin(pointBounds.left, pointBounds.right, 0, this.bounds.width),
				pair[1].y.linlin(pointBounds.top, pointBounds.bottom, 0, this.bounds.height)
			);
			Pen.line(scaledPoint1, scaledPoint2);
		};
		Pen.stroke;
	}

	thresh2 {
		arg x;
		^x.abs.thresh(this.errorAllowance) * x.sign;
	}

	giftwrap {
		arg points;
		var leftmost, perimeter, p_curr, p_cand, prev_theta;
		// var points_temp = List[];

		/*// remove duplicates (within error)
		for(0,points.size-1) {
			|i|
			var flag = false;
			if(i != (points.size - 1)) {
				for(i+1, points.size-1) {
					|j|
					var pa = points[i], pb = points[j];
					flag = flag || (pa.x.equalWithPrecision(pb.x,errorAllowance) &&
						pa.y.equalWithPrecision(pb.y,errorAllowance));
				};
			};
			flag.not.if {points_temp = points_temp.add(points[i])};
		};
		points = points_temp;*/

		leftmost = points.sort({|a,b| (a.x==b.x).if {a.y>b.y} {a.x<b.x}}).first;
		perimeter = [leftmost];
		p_cand = nil;
		prev_theta = this.giftwrap_theta(p_curr, nil);

		while {(p_cand != perimeter.first) || (perimeter.size < 2)} {
			p_cand = points[0];
			for(1, points.size-1) {
				|j|
				if((p_cand == perimeter.last) || this.giftwrap_leftOfLine(perimeter.last, p_cand, points[j], prev_theta)) {
					p_cand = points[j];
				}
			};
			prev_theta = this.giftwrap_theta(p_cand, perimeter.last);
			perimeter = perimeter.add(p_cand);
		};

		perimeter = perimeter.slide(2,1).clump(2);

		^perimeter;
	}

	giftwrap_theta {
		arg p_curr, p_prev;
		^p_prev.isNil.if {
			0.5pi;
		} {
			(p_curr - p_prev).theta;
		}
	}

	giftwrap_leftOfLine {
		arg p_curr, p_cand, p_test, prev_theta;

		var theta_cand = prev_theta - (p_cand - p_curr).theta;
		var theta_test = prev_theta - (p_test - p_curr).theta;

		// "cand: %\ttheta_cand: %".format(p_cand, theta_cand).postln;
		// "test: %\ttheta_test: %".format(p_test, theta_test).postln;

		^((theta_test % 2pi) < (theta_cand % 2pi)) && (p_test != p_curr);
	}





}