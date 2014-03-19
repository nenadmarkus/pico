//
var n3_find_faces = Module.cwrap('find_faces', 'number', (['number', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number']));

// allocate memory
var ppixels = Module._malloc(640*480);
var pixels = new Uint8Array(Module.HEAPU8.buffer, ppixels, 640*480);

var maxndetections = 1024;

prs = Module._malloc(4*maxndetections)
var rs = new Float32Array(Module.HEAPU8.buffer, prs, maxndetections);

pcs = Module._malloc(4*maxndetections)
var cs = new Float32Array(Module.HEAPU8.buffer, pcs, maxndetections);

pss = Module._malloc(4*maxndetections)
var ss = new Float32Array(Module.HEAPU8.buffer, pss, maxndetections);

pqs = Module._malloc(4*maxndetections)
var qs = new Float32Array(Module.HEAPU8.buffer, pqs, maxndetections);

//
function process_frame(mCanvas)
{
	var mContext = mCanvas.getContext("2d"),
		mNewData = mContext.getImageData(0, 0, mCanvas.width, mCanvas.height),
		rgba = mNewData.data
	;

	for(i=0; i<rgba.length; i+=4)
	{
		average = 0.299*rgba[i] + 0.587*rgba[i+1] + 0.114*rgba[i+2];
		pixels[i/4] = average;
	}

	//
	var ndetections = n3_find_faces(prs, pcs, pss, pqs, maxndetections, ppixels, mCanvas.height, mCanvas.width, mCanvas.width);

	//console.log(ndetections);

	// draw results
	var centerX, centerY, radius;
	for(i=0; i<ndetections; ++i)
	{
		if(qs[i]>3.0)
		{
			centerX = cs[i];
			centerY = rs[i];
			radius = ss[i]/2;

			mContext.beginPath();
			mContext.arc(centerX, centerY, radius, 0, 2*Math.PI, false);
			mContext.lineWidth = 3;
			mContext.strokeStyle = 'red';
			mContext.stroke();
		}
	}
}

var video = document.createElement('video');

function startStream(stream)
{
	video.addEventListener('canplay', function DoStuff()
	{
		video.removeEventListener('canplay', DoStuff, true);
		setTimeout(function() {
			video.play();

			canvas.width = video.videoWidth;
			canvas.height = video.videoHeight;

			drawToCanvas();
		}, 3000);
	}, true);

	var domURL = window.URL || window.webkitURL;
	video.src = domURL ? domURL.createObjectURL(stream) : stream;

	video.play();
}
function deniedStream()
{
	alert("Camera access denied!");
}
function errorStream(e)
{
	if (e)
	{
		alert(e);
	}
}
function drawToCanvas()
{
	window.requestAnimationFrame(drawToCanvas);

	canvas.getContext('2d').drawImage(video, 0, 0, canvas.width, canvas.height);

	process_frame(canvas);
}

// initialize video input
var i=0, vendors=['ms', 'moz', 'webkit', 'o'];

while (i<vendors.length && !window.requestAnimationFrame)
{
	window.requestAnimationFrame = window[vendors[i] + 'RequestAnimationFrame'];
	i++;
}

if(!window.requestAnimationFrame)
{
	alert("RequestAnimationFrame mechanism is not supported by this browser.");
}

//
navigator.getUserMedia_ = navigator.getUserMedia || navigator.webkitGetUserMedia || navigator.mozGetUserMedia || navigator.msGetUserMedia;

try
{
	navigator.getUserMedia_({
		video: true,
		audio: false
	}, startStream, deniedStream);
}
catch (e)
{
	try
	{
		navigator.getUserMedia_('video', startStream, deniedStream);
	}
	catch (e)
	{
		errorStream(e);
	}
}

video.loop = video.muted = true;
video.load();
