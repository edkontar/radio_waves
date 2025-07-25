function validateFrequencyInput(input) {
  // Parse the input value as a floating-point number
  const enteredValue = parseFloat(input.value);

  // Check if the entered value is less than the minimum allowed
  if (enteredValue < parseFloat(input.min)) {
    input.setCustomValidity('Value must be greater than or equal to ' + input.min);
  } else {
    input.setCustomValidity('');
  }
}

// JavaScript function to toggle checkboxes
function toggleCheckbox(checkbox) {
  var checkboxes = document.getElementsByName(checkbox.name);
  checkboxes.forEach(function (item) {
    if (item !== checkbox) {
      item.checked = false;
    }
  });
}

// Function to get the appropriate array based on the checkbox status
async function getArraysBasedOnCheckbox() {
  
  const useFundamental = document.getElementById('checkbox1').checked;
  const useHarmonic = document.getElementById('checkbox2').checked;

  // Determine the JSON file based on checkbox status
  let jsonFile;
  if (useFundamental && !useHarmonic) {
    jsonFile = 'https://raw.githubusercontent.com/edkontar/radio_waves/main/results/new_all_results_F_final.json';
  } else if (!useFundamental && useHarmonic) {
    jsonFile = 'https://raw.githubusercontent.com/edkontar/radio_waves/main/results/new_all_results_H_final.json';
  } else {
    // Handle the case where both checkboxes are checked or none are checked
    console.error('Invalid checkbox combination');
    return {}; // Return some default values or handle the error as needed
  }

  try {
    // Fetch the JSON data
    const response = await fetch(jsonFile);
    
    if (!response.ok) {
      console.error('Error fetching JSON:', response.status, response.statusText);
      return {}; // Return some default values or handle the error as needed
    }

    // Access the arrays directly from the loaded JSON data
    const json = await response.json();

    return {
      fmhz: json.FMHZ,
      eps: json.EPS,
      anis: json.ANIS,
      r_init: json.R_INIT,
      r_shift: json.R_SHIFT,
      s_fwhm: json.S_FWHM,
      t_1e_decay: json.T_1E_AFTER,
    };
  } catch (error) {
    console.error('Error fetching JSON:', error);
    return {}; // Return some default values or handle the error as needed
  }
}


// Generate 1D arrays
function range(start, end, numElements) {
  const nums = [];
  const step = (end - start) / (numElements - 1);
  for (let i = 0; i < numElements; i++) {
    nums.push(start + (step * i));
  }
  return nums;
};

function interpolate(x, y, xValue) {
  // Check if xValue is less than or equal to the first value in x
  if (xValue <= x[0]) {
    return y[0];
  }
  // Check if xValue is greater than or equal to the last value in x
  if (xValue >= x[x.length - 1]) {
    return y[y.length - 1];
  }
  // Find the index i where x[i] <= xValue < x[i+1]
  let i = 0;
  while (x[i] < xValue) {
    i++;
  }
  // Perform linear interpolation
  const x0 = x[i - 1];
  const x1 = x[i];
  const y0 = y[i - 1];
  const y1 = y[i];
  // Calculate the interpolated value
  const interpolatedY = y0 + ((xValue - x0) * (y1 - y0)) / (x1 - x0);
  return interpolatedY;
};

// Generate 2D Gaussian 
function makeGaussian(amplitude, x0, y0, sigmaX, sigmaY) {
  return function(x, y) {
    var exponent = -(
        Math.pow(x - x0, 2) / (2 * Math.pow(sigmaX, 2)) +
        Math.pow(y - y0, 2) / (2 * Math.pow(sigmaY, 2))
      );
    return amplitude * Math.pow(Math.E, exponent);
  }
};

// find all indexes of a given value in an array
function getAllIndexes(arr, val) {
  var indexes = [], i;
  for(i = 0; i < arr.length; i++)
      if (arr[i] === val)
          indexes.push(i);
  return indexes;
};

// rotate 1D array
function rotate1D(arr, angleInRadians) {
  return arr.map(function (value) {
    return value * Math.cos(angleInRadians);
  });
}

// rotate 1D x array around yaxis
function rotateX_Y(arr, angleInRadians, shift) {
  return arr.map(function (value) {
    return value * Math.cos(angleInRadians) + shift * Math.sin(angleInRadians);
  });
}

// find all data matching a given eps and anis combination
function fMatchUserEpsAnis(fmhz0, eps0, anis0, chosenEps, chosenAnis) {
  const indices = [];
  // Iterate through eps and anis arrays to find matching indices
  for (let i = 0; i < eps0.length; i++) {
    if (eps0[i] === chosenEps && anis0[i] === chosenAnis) {
      indices.push(i);
    }
  }
  // Use the indices to get the corresponding fmhz values
  const correspondingFrequencies = indices.map(index => fmhz0[index]);
  return { indices, correspondingFrequencies };
};

// Function to update the frequency slider value display
function updateFreqSliderValue(sliderId, displayId) {
  const slider = document.getElementById(sliderId);
  const display = document.getElementById(displayId);
  display.textContent = slider.value + " MHz";
}

// Function to update the helio angle slider value display
function updateAngSliderValue(sliderId, displayId) {
  const slider = document.getElementById(sliderId);
  const display = document.getElementById(displayId);
  display.textContent = slider.value + " deg";
}

// Function to tile 1D array to 2D array
function tileVertically(originalArray, numRows, numCols) {
  if (originalArray.length !== numRows) {
    throw new Error("Invalid dimensions for tiling.");
  }

  const resultArray = [];

  for (let i = 0; i < numRows; i++) {
    const row = Array.from({ length: numCols }, () => originalArray[i]);
    resultArray.push(row);
  }

  return resultArray;
}

function degToRad(deg) {
  return deg * (Math.PI / 180);
}

// for creating helio grid lat,lon to arcsec from solar center
function convertToXYZ(lat, lon, radius, xRotationDegrees) {
  const latRad = degToRad(lat);
  const lonRad = degToRad(lon);
  const xRotationRad = degToRad(xRotationDegrees);

  const y = radius * Math.sin(latRad);
  const x = radius * Math.cos(latRad) * Math.sin(lonRad);
  const z = -radius * Math.cos(latRad) * Math.cos(lonRad);

  // Rotation matrix around x-axis
  const rotatedY = y * Math.cos(xRotationRad) - z * Math.sin(xRotationRad);
  const rotatedZ = y * Math.sin(xRotationRad) + z * Math.cos(xRotationRad);

  return { x, y: rotatedY, z: rotatedZ };
}

// Run the updateSliderValue function when the page loads
window.addEventListener('load', function () {
  //updateFreqSliderValue("frequency", "frequencyValue");
  updateAngSliderValue("helioAngle", "helioAngleValue");
});

// Add an event listener to the slider to update the value display in real-time
//document.getElementById("frequency").addEventListener('input', function () {
//  updateFreqSliderValue("frequency", "frequencyValue");
//});
document.getElementById("helioAngle").addEventListener('input', function () {
  updateAngSliderValue("helioAngle", "helioAngleValue");
});

let debounceTimeout;

function plotGraphs() {

  // Clear any existing timeout to prevent previous inputs from triggering the function
  clearTimeout(debounceTimeout);

  // Set a new timeout to run the function after 1000ms (1 second)
  debounceTimeout = setTimeout(function () {

    const frequencyInput = document.getElementById("frequency");
    // Parse the entered frequency value
    const f_user = parseFloat(frequencyInput.value);

    // Validate the entered frequency
    if (isNaN(f_user) || f_user < parseFloat(frequencyInput.min)) {
      // Display an alert or handle the invalid input in some way
      alert('Frequency must be greater than or equal to ' + frequencyInput.min);
      
      // Reset the value to the minimum allowed
      frequencyInput.value = frequencyInput.min;

      return; // Exit the function early
    }

    // Check if the entered frequency exceeds the maximum allowed
    if (f_user > parseFloat(frequencyInput.max)) {
      // Set the value back to the maximum allowed
      frequencyInput.value = frequencyInput.max;
      
      // Optionally, you can display an alert or handle the situation in some way
      alert('Frequency must be between ' + frequencyInput.min + ' and ' + frequencyInput.max);
      
      return; // Exit the function early
    }

    var eps_user = parseFloat(document.getElementById("epsilon").value);
    const anis_user = parseFloat(document.getElementById("anisotropy").value);
    const helioAng_user = parseFloat(document.getElementById("helioAngle").value);
    const helioAngRad_user = (helioAng_user * Math.PI)/180;
  
    // match rounded user entry to data values
    if (eps_user === 0.7) {eps_user = 0.709};
    if (eps_user === 1.4) {eps_user = 1.41};

    async function run() {

        const {
          r_init,
          r_shift,
          fmhz,
          eps,
          anis,
          s_fwhm,
          t_1e_decay,
        } = await getArraysBasedOnCheckbox();

      // User value index search
      const fidx = fmhz.findIndex(num => num > f_user)
      const fidxs = getAllIndexes(fmhz,fmhz[fidx])
      // Initialize a variable to store the final index
      var finalIndex = -1;
      // Loop through the matching indices
      for (var j = 0; j < fidxs.length; j++) {
        var index = fidxs[j];
        // Check if eps and anis have the desired values at this index
        if (eps[index] === eps_user && anis[index] === anis_user) {
          finalIndex = index;
        break; // Found the desired index, so exit the loop
        }
      };

      // Get data matching user chosen eps and anis
      const freqsObj = fMatchUserEpsAnis(fmhz, eps, anis, eps_user, anis_user);
      const freqs = freqsObj.correspondingFrequencies;
      const sizesDeg_FWHM = freqsObj.indices.map(index => s_fwhm[index]);
      const sizesDeg = sizesDeg_FWHM.map(value => value/2.355); // converts from FWHM->stdev
      const sizes = sizesDeg.map(value => value*3600); // converts from deg->arcsec
      const decay = freqsObj.indices.map(index => t_1e_decay[index]);
      const rr = freqsObj.indices.map(index => r_init[index]);
      const rshift = freqsObj.indices.map(index => r_shift[index]);
      
      // Interpolate at user frequency
      var sint = interpolate(freqs,sizes,f_user); // arcsec
      var sint_fwhm = sint/3600*2.355; // for size,freq plot which is always in deg
      const dint = interpolate(freqs,decay,f_user); // s
      const rint = interpolate(freqs,rr,f_user); // Rsun
      const shiftint = interpolate(freqs,rshift,f_user); // Rsun

      var xlabel = 'X [arcsec]';
      var ylabel = 'Y [arcsec]';
      var xRange = 6 * sint;
      var yRange = xRange;
      var slabel = 'arcsec';
      var rsun = 950; // arcsec
      
      if (xRange > 10000) {
        xRange = xRange/3600;
        yRange = xRange;
        xlabel = 'X [Degrees]';
        ylabel = 'Y [Degrees]';
        sint = sint/3600; // arcsec->deg
        slabel = 'degrees';
        rsun = rsun/3600; // arcsec->deg
      };

      // Generate plot data
      const gaussian = makeGaussian(1, 0, 0, sint, sint);

      // create helio grid
      const grid_spacing = 10;
      const gv = range(0, 36, 37).map(value => value * grid_spacing - 90).flat();
      const np = gv.length;
      const lon = tileVertically(gv, np, np);
      const lat = lon[0].map((val, index) => lon.map(row => row[index]).reverse()); // rotate lon by 90 deg      
      const latLonCoordinates = [];
  
      for (let i = 0; i < np; i++) {
        for (let j = 0; j < np; j++) {
          const latVal = lat[i][j];
          const lonVal = lon[i][j];
          const { x, y, z } = convertToXYZ(latVal, lonVal, rsun, 0);
          latLonCoordinates.push({ lat: latVal, lon: lonVal, x, y, z });
        }
      }

      // create x,y arrays for image
      const numXYPoints = 300;
      var x = range(-xRange, xRange, numXYPoints);
      var y = x;
      var z = shiftint*rsun;
      var xc = 0;
      var im0 = [];
      for (let i = 0; i < x.length; i++) {
        im0.push(y.map(y => gaussian(x[i], y))); 
      };

      // Rotate around y-axis
      x = rotateX_Y(x, helioAngRad_user, z);
      xc = shiftint*rsun*Math.sin(helioAngRad_user);

      // Find the minimum and maximum values in the data
      // Iterates because Safari has issues with max stack using Math.min()
      const flattenedZ0min = im0.flat();
      let minZ = flattenedZ0min[0];
      for (let i = 1; i < flattenedZ0min.length; i++) {
        if (flattenedZ0min[i] < minZ) {
          minZ = flattenedZ0min[i];
        }
      };

      const flattenedZ0max = im0.flat();
      let maxZ = flattenedZ0max[0];
      for (let i = 1; i < flattenedZ0max.length; i++) {
        if (flattenedZ0max[i] > maxZ) {
          maxZ = flattenedZ0max[i];
        }
      };

      const minS = Math.min(...sizesDeg); // deg
      const maxS = Math.max(...sizesDeg); // deg

      // Normalize the data to the range [0, 1]
      const im = im0.map(row => row.map(value => (value - minZ) / (maxZ - minZ)));

      // Calculate the FWHM value
      // interate to find max value in data because Safari fails using Math.max()
      const flattenedZ = im.flat();
      let maxAmplitude = flattenedZ[0];
      for (let i = 1; i < flattenedZ.length; i++) {
        if (flattenedZ[i] > maxAmplitude) {
          maxAmplitude = flattenedZ[i];
        }
      }
      //const maxAmplitude = Math.max(...z.flat()); // Find the maximum value in the data
      const fwhmThreshold50 = maxAmplitude * 0.50;
      const fwhmThreshold70 = maxAmplitude * 0.70;
      const fwhmThreshold90 = maxAmplitude * 0.90;

      // Create a range of contour levels from 0 to 1
      const contourLevels = Array.from({ length: 11 }, (_, i) => i / 10); // 11 levels from 0 to 1
      
      // Create a contour trace with the defined contour levels
      const contourTrace = {
        type: 'heatmap',
        x: x,
        y: y,
        z: im,
        colorscale: 'YlGnBu',
        colorbar: {
          title: 'Normalised Intensity',
          titleside: 'right',
          len: 1.06,
        },
        hoverinfo: 'skip',
        zmin: 0, // Set the minimum value for the color scale
        zmax: 1, // Set the maximum value for the color scale
      };

      // Create a contour trace at the 50% level
      const fwhmTraces = [
        {
          type: 'contour',
          x: x,
          y: y,
          z: im,
          showlegend: false,
          hoverinfo: 'skip',
          name: 'FWHM',
          contours: {
            start: fwhmThreshold50.toFixed(1),
            end: fwhmThreshold50.toFixed(1),
            size: 0, // Set to 0 to only plot the 50% contour line
            coloring: 'none',
            showlabels: true,
            labelfont: {
              color: 'black'
            },
            showlines: true,
            line: {
              dash: 'dash', // Make the contour line dashed
              width: 4, // Adjust the line width as needed
              color: 'black', // Set the line color to black
            },
          },
        },
        {
          type: 'contour',
          x: x,
          y: y,
          z: im,
          showlegend: false,
          hoverinfo: 'skip',
          name: 'FWHM',
          contours: {
            start: fwhmThreshold70.toFixed(1),
            end: fwhmThreshold70.toFixed(1),
            size: 0, // Set to 0 to only plot the 50% contour line
            coloring: 'none',
            showlabels: true,
            labelfont: {
              color: 'black'
            },
            showlines: true,
            line: {
              dash: 'dash', // Make the contour line dashed
              width: 4, // Adjust the line width as needed
              color: 'black', // Set the line color to black
            },
          },
        },
        {
          type: 'contour',
          x: x,
          y: y,
          z: im,
          showlegend: false,
          hoverinfo: 'skip',
          name: 'FWHM',
          contours: {
            start: fwhmThreshold90.toFixed(1),
            end: fwhmThreshold90.toFixed(1),
            size: 0, // Set to 0 to only plot the 50% contour line
            coloring: 'none',
            showlabels: true,
            labelfont: {
              color: 'black'
            },
            showlines: true,
            line: {
              dash: 'dash', // Make the contour line dashed
              width: 4, // Adjust the line width as needed
              color: 'black', // Set the line color to black
            },
          },
        },
      ];

      const solarGridTrace = {
        type: 'scatter',
        mode: 'lines',
        x: [],
        y: [],
        line: {
          //dash: 'dot',
          opacity: 0.1,
          color: 'white',
          width: 0.025,
        },
        hoverinfo: 'skip',
        name: ' ',
      };

      // Create a 2D grid-like pattern by connecting points in the correct order
      for (let i = 0; i < np; i++) {
        const xSubset = latLonCoordinates.slice(i * np, (i + 1) * np).map(coord => coord.x);
        const ySubset = latLonCoordinates.slice(i * np, (i + 1) * np).map(coord => coord.y);

        // Set the style for the additional lines
        const additionalLines = {
          type: 'scatter',
          mode: 'lines',
          x: xSubset.concat(xSubset[0]),
          y: ySubset.concat(ySubset[0]),
          line: {
            //dash: 'dot',
            opacity: 0.1,
            color: 'white',
            width: 0.025,
          },
          hoverinfo: 'skip',
          name: ' ',
        };

        // Concatenate the additional lines to the trace
        solarGridTrace.x = solarGridTrace.x.concat(additionalLines.x); 
        solarGridTrace.y = solarGridTrace.y.concat(additionalLines.y);
      }

      // Create additional lines along the other axis
      for (let j = 0; j < np; j++) {
        const xSubset = latLonCoordinates.filter((_, index) => index % np === j).map(coord => coord.x);
        const ySubset = latLonCoordinates.filter((_, index) => index % np === j).map(coord => coord.y);

        // Set the style for the additional lines
        const additionalLines = {
          type: 'scatter',
          mode: 'lines',
          x: xSubset.concat(xSubset[0]),
          y: ySubset.concat(ySubset[0]),
          line: {
            //dash: 'dot',
            opacity: 0.1,
            color: 'white',
            width: 0.025,
          },
          hoverinfo: 'skip',
          name: ' '
        };

        // Concatenate the additional lines to the trace
        solarGridTrace.x = solarGridTrace.x.concat(additionalLines.x);
        solarGridTrace.y = solarGridTrace.y.concat(additionalLines.y);
      }

      var data = [contourTrace, ...fwhmTraces, solarGridTrace];

      // Plot Sun
      const numPoints = 1000; // Generate an array of angles from 0 to 2pi
      const circle_angles = Array.from({ length: numPoints }, (_, i) => (i / (numPoints - 1)) * 2 * Math.PI);
      // Calculate x and y coordinates of the points on the circle
      const sunX = circle_angles.map(angle => rsun * Math.cos(angle));
      const sunY = circle_angles.map(angle => rsun * Math.sin(angle));
      // Create the Plotly trace for the circle
      const sunTrace = {
        x: sunX,
        y: sunY,
        mode: 'lines',
        line: {
          color: 'white', // Adjust the color as needed
          width: 2,     // Adjust the line width as needed
        },
        hoverinfo: 'skip',
        name: ' ', // Legend label for the circle
      };

      // Create size against freq plot
      const sizeFreqTrace = {
        x: freqs,
        y: sizesDeg_FWHM,
        type: 'scatter',
        mode: 'lines+markers',
        showlegend: false
      };

      const intSizeYOverlayTrace = {
        x: [f_user,f_user],
        y: [0,sint_fwhm],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      const intSizeXOverlayTrace = {
        x: [0,f_user],
        y: [sint_fwhm,sint_fwhm],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      data2 = [sizeFreqTrace, intSizeXOverlayTrace, intSizeYOverlayTrace];

      // Create position plot
      const xpos = rint*Math.sin(helioAngRad_user);
      const zpos = rint*Math.cos(helioAngRad_user);
      const posTrace = {
        x: [xpos],
        y: [zpos],
        type: 'scatter',
        mode: 'markers',
        marker: {
          symbol: 'x',
          size: 8,  // smaller size makes lines thinner visually
          color: 'red',
        },
        showlegend: false,
        name: 'Source'
      };

      // Sun for pos plot
      const sunPosTrace = {
        x: circle_angles.map(angle => Math.cos(angle)),
        y: circle_angles.map(angle => Math.sin(angle)),
        type: 'scatter',
        mode: 'lines+markers',
        marker: {
          symbol: 'circle',
          color: '#FFD700',
          size: 3
        },
        showlegend: true,
        name: 'Sun',
        hoverinfo: 'skip'
      };

      data3 = [sunPosTrace, posTrace];

      // Trace for radii circles in position plot
      const radii = [10, 20, 30, 40, 50, 60, 70, 80];
      const labelAngle = -45; // Desired angle for the labels (in degrees)
      const circleTraces = radii.map(radius => ({
        x: circle_angles.map(angle => radius * Math.cos(angle)),
        y: circle_angles.map(angle => radius * Math.sin(angle)),
        mode: 'lines',
        line: {
            dash: 'dash',
            color: 'grey',
            width: 1,
        },
        hoverinfo: 'skip',
        showlegend: false,
        name: ' ', // Legend label for the circle
      }));
      
      // Create a separate trace for the text labels
      const circleLabelAnnotations = radii.map((radius, index) => {
        const labelX = radius * Math.cos(labelAngle * (Math.PI / 180)); // Convert angle to radians
        const labelY = radius * Math.sin(labelAngle * (Math.PI / 180)); // Convert angle to radians
        return {
          x: labelX,
          y: labelY,
          text: `${radius}`,
          xref: 'x',
          yref: 'y',
          showarrow: false,
          font: {
            color: 'black',
            size: 10,
          },
          bgcolor: 'rgba(255, 255, 255, 0.5)', // Transparent white background for the label
        };
      });

      // Trace for last scattering radius in position plot
      const lastScatTrace = {
        x: circle_angles.map(angle => shiftint * Math.cos(angle)),
        y: circle_angles.map(angle => shiftint * Math.sin(angle)),
        mode: 'lines',
        line: {
          dash: 'solid',
          color: 'red',
          width: 1,
        },
        hoverinfo: 'skip',
        showlegend: false,
        name: ' ', // Legend label for the circle
      };
      
      // Combine all the circle traces into a single data array
      const data3b = [...circleTraces, lastScatTrace];

      // Create size against freq plot
      const decayTrace = {
        x: freqs,
        y: decay,
        type: 'scatter',
        mode: 'lines+markers',
        showlegend: false
      };

      const intDecayYOverlayTrace = {
        x: [f_user,f_user],
        y: [0,dint],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      const intDecayXOverlayTrace = {
        x: [0,f_user],
        y: [dint,dint],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      data4 = [decayTrace, intDecayXOverlayTrace, intDecayYOverlayTrace];

      // Create distance against freq plot
      const distTrace = {
        x: freqs,
        y: rshift.map(val => val - 1),
        type: 'scatter',
        mode: 'lines+markers',
        showlegend: false
      };

      const intDistYOverlayTrace = {
        x: [f_user,f_user],
        y: [0,shiftint-1],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      const intDistXOverlayTrace = {
        x: [0,f_user],
        y: [shiftint-1,shiftint-1],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      const fpeDistXOverlayTrace = {
        x: [0,f_user],
        y: [rint-1,rint-1],
        mode: 'lines',
        line: {color: 'black', width: 1, dash: 'dash'},
        hoverinfo: 'skip',
        name: ' ',
        showlegend: false
      };

      const fpeTrace = {
        x: freqs,
        y: rr.map(val => val - 1),
        type: 'scatter',
        mode: 'lines', // or 'lines+markers' if you want markers too
        name: 'rint', // label for legend
        showlegend: false,
        line: { width: 2, color: 'black' } // optional styling
      };

      data5 = [distTrace, intDistXOverlayTrace, intDistYOverlayTrace, fpeTrace, fpeDistXOverlayTrace];

      // Set image plot style
      var layout1 = {
        width: 500,  // set width
        height: 485, // set height
        aspectratio: { x: 1, y: 1 }, // Set equal aspect ratio for x and y
        plot_bgcolor: '#800026',
        title: {
          text: 'Apparent Source Image',
          y: 0.85
        },
        xaxis: {
          title: xlabel,
          range: [-xRange/2+xc, xRange/2+xc],
          showgrid: false,
          zeroline: false,
        },
        yaxis: {
          title: {
            text: ylabel,
            standoff: 20
          },
          range: [-yRange/2, yRange/2],
          showgrid: false,
          zeroline: false,
        },
        showlegend: true,
        legend: {
          //x: -3,
          y: -20,
        },
        annotations: [
        {
          x: -0.1,
          y: 1.3,
          xref: 'paper',
          yref: 'paper',
          text: 'Source Heliocentric Distance: '+rint.toFixed(2)+' R<sub>☉</sub>',
          showarrow: false,
        },
        {
          x: -0.1,
          y: 1.25,
          xref: 'paper',
          yref: 'paper',
          text: 'Frequency: '+f_user.toFixed(1)+' MHz',
          showarrow: false,
        },
        //{
        // x: -0.1,
        //  y: 1.2,
        //  xref: 'paper',
        //  yref: 'paper',
        //  text: 'FWHM Major Size: '+(sint*2.355).toFixed(1)+' '+slabel,
        //  showarrow: false,
        //},
        ],
        dragmode: 'pan', // disables zooming
      };

      // Set sizeFreq plot style
      var layout2 = {
        width: 500,  // set width
        height: 500, // set height
        aspectratio: { x: 1, y: 1 }, // Set equal aspect ratio for x and y
        xaxis: {
          title: 'Frequency [MHz]',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          type: 'log',
          linewidth: 2,
          range: [-1, 2.53],
          showline: true,
        },
        yaxis: {
          title: 'FWHM Major Size [degrees]',
          type: 'log',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          linewidth: 2,
          autorange: true,
          showline: true, 
        },
        dragmode: false,
        // Add annotation for FWHM Major Size value
        annotations: [
          {
            x: 1,
            y: 1,
            xref: 'paper',
            yref: 'paper',
            text: 'FWHM Major Size: ' + (sint * 2.355).toFixed(1) + ' ' + slabel,
            showarrow: false,
          },
        ],
      };
      
      // last scattering radius annotation
      const additionalAnnotations = [
        {
          x: 0,
          y: 1,
          xref: 'paper',
          yref: 'paper',
          text: 'Radius of Last Scattering',
          showarrow: false,
          showlegend: false,
          font: {
            color: 'red', // Set the color property to 'red'
          },
        },
      ];
      // combine annotations for position plot
      const allAnnotations = [...circleLabelAnnotations, ...additionalAnnotations];

      // Set position plot style
      var layout3 = {
        width: 480,  // set width
        height: 500, // set height
        aspectratio: { x: 1, y: 1 }, // Set equal aspect ratio for x and y
        showlegend: true,
        annotations: allAnnotations,
        pad: 0,
        xaxis: {
          title: 'X [R<sub>☉</sub>]',
          ticks: 'inside',
          //tickvals: [-80, -60, -40, -20, 0, 20, 40, 60, 80],
          dtick: Math.round((shiftint * 1.2 * 2) / 5),
          linewidth: 2,
          range: [-shiftint*1.2, shiftint*1.2],
          showline: true,
        },
        yaxis: {
          title: 'Z [R<sub>☉</sub>]',
          ticks: 'inside',
          //tickvals: [-80 -60, -40, -20, 0, 20, 40, 60, 80],
          dtick: Math.round((shiftint * 1.2 * 2) / 5),
          linewidth: 2,
          range: [shiftint*1.2, -shiftint*1.2],
          showline: true, 
        },
        dragmode: false,
      };

      // Set decay plot style
      var layout4 = {
        width: 500,  // set width
        height: 500, // set height
        aspectratio: { x: 1, y: 1 }, // Set equal aspect ratio for x and y
        xaxis: {
          title: 'Frequency [MHz]',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          type: 'log',
          linewidth: 2,
          range: [-1, 2.53],
          showline: true,
        },
        yaxis: {
          title: '1/e Decay Time [s]',
          type: 'log',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          linewidth: 2,
          autorange: true,
          showline: true, 
        },
        dragmode: false,
        // Add annotation for decay value
        annotations: [
          {
            x: 1,
            y: 1,
            xref: 'paper',
            yref: 'paper',
            text: '1/e Decay: ' + (dint).toFixed(2) + ' s',
            showarrow: false,
          },
        ],
      };

      // Set distance plot style
      var layout5 = {
        width: 500,  // set width
        height: 500, // set height
        aspectratio: { x: 1, y: 1 }, // Set equal aspect ratio for x and y
        xaxis: {
          title: 'Frequency [MHz]',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          type: 'log',
          linewidth: 2,
          range: [-1, 2.53],
          showline: true,
        },
        yaxis: {
          title: 'Distance [r/R<sub>☉</sub> - 1]',
          type: 'log',
          dtick: 1, // show only powers of 10 labels
          ticks: 'inside',
          linewidth: 2,
          autorange: true,
          showline: true, 
        },
        dragmode: false,
        // Add annotation for decay value
        annotations: [
          {
            x: 1,
            y: 0.95,
            xref: 'paper',
            yref: 'paper',
            text: 'Apparent Source Height: ' + (shiftint-1).toFixed(2) + ' R<sub>☉</sub>',
            showarrow: false,
          },
          {
            x: 1,
            y: 1,
            xref: 'paper',
            yref: 'paper',
            text: 'True Source Height: ' + (rint-1).toFixed(2) + ' R<sub>☉</sub>',
            showarrow: false,
          },
          {
            x: 1,
            y: 0.9,
            xref: 'paper',
            yref: 'paper',
            text: 'Shift: ' + (shiftint-rint).toFixed(2) + ' R<sub>☉</sub>',
            showarrow: false,
          },
        ],
      };

      const config = {
        displayModeBar: true,
        toImageButtonOptions: {
          format: 'png',
          filename: 'scattered_image_'+f_user.toPrecision(3)+'MHz_'+helioAng_user.toPrecision(2)+'deg_eps'+eps_user.toPrecision(2)+'_anis'+anis_user.toPrecision(2),
          scale: 2
        },
        modeBarButtonsToRemove: ['autoScale2d','toggleSpikelines','hoverClosestCartesian','hoverCompareCartesian']
      };

      const config2 = {
        displayModeBar: false,
      };

      // Draw plots
      Plotly.newPlot('imageDiv', data, layout1, config);
      Plotly.addTraces('imageDiv', sunTrace);
      Plotly.newPlot('sizeFreqDiv', data2, layout2, config2);
      Plotly.newPlot('posDiv', data3, layout3, config2);
      Plotly.addTraces('posDiv', data3b);
      Plotly.newPlot('decayDiv', data4, layout4, config2);
      Plotly.newPlot('distDiv', data5, layout5, config2);
    }
           
  // Call the async function
  run();

  }, 1000); // Wait for 1000ms (1 second) before executing the code
}


// Run the plotGraphs function when the page loads
window.addEventListener('load', plotGraphs);

// Run the plotGraphs function when any form input value changes
const formInputs = document.querySelectorAll('form input, form select');
formInputs.forEach((input) => {
  input.addEventListener('input', plotGraphs);
  // Add an event listener for the "keydown" event to prevent form submission on Enter key press
  input.addEventListener('keydown', function (event) {
    if (event.key === 'Enter') {
      event.preventDefault();
      plotGraphs();
    }
  });
});
