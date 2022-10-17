function findzero(arr) {
    var id = [0] //becomes array, just number w/o []
    var k = 0
    for (let i = 0; i < (arr.length - 1); i++) {
        if (Math.sign(arr[i]) == 0) {
            id[k] = i
            k++
        } else if (Math.sign(arr[i]) != Math.sign(arr[i + 1])) {
            id[k] = i
            k++
        }
    }
    return id
}

function linspace(min, max, nel) {
    nel = nel - 1
    var xMax = max; //max x   
    var xMin = min; //min x from previous section
    var x = [...Array(nel + 1).keys()];
    x = x.map(a => a * ((xMax - xMin) / nel) + xMin)
    return x
}

function getWidth(data) {
    const { h, x_p, y_p } = data;
    y_c = linspace(0, h, h + 1)
    const width = []
    for (let j = 0; j < y_c.length; j++) {
        let xCrossing = [];
        let w = [];
        //loop over relevant segments to find intersect
        for (let i = 0; i < x_p.length - 1; i++) {
            //check if segment cross the present yLine
            if ((((y_p[i] > y_c[j]) + (y_p[i + 1] > y_c[j])) & 1) > 0.5) {
                //find and store x-value of crossing
                if (Math.abs(x_p[i] - x_p[i + 1]) < 0.0001) {
                    xIntersect = x_p[i];
                } else {
                    const k = (y_p[i] - y_p[i + 1]) / (x_p[i] - x_p[i + 1]);
                    const m = y_p[i] - k * x_p[i]
                    xIntersect = (y_c[j] - m) / k
                }
                xCrossing = [...xCrossing, xIntersect]
                //console.log(xCrossing)

            }

        }
        for (let k = 0; k < xCrossing.length; k = +2) {
            w = [...w, Math.abs(xCrossing[k] - xCrossing[k + 1])]
        }
        width[j] = w.reduce((a, b) => a + b, 0);
    }
    //ax = document.getElementById('myAx3')
    //Plotly.newPlot(ax, [{ x: x_p, y: y_p }, { x: width, y: y_c }])
    return [y_c, width]
}

function calc() {
    data = {}
    let groups = $('.geom').each(function () {
        let id = $(this).attr('id');
        // create field name from id and index. e.g. fld1_1
        const field = `${id}`;
        data = {
            ...data,  // spread current items
            [field]: parseInt($(this).val())  // set new input value
        }

    })
    data = { ...data, b: Math.max(data.b_f1, data.b_f2) }
    const { h, t_f2, t_f1, t_w, b_f1, b_f2, b } = data

    //define cross section polygon
    x_p = [b / 2 - b_f1 / 2, b / 2 + b_f1 / 2, b / 2 + b_f1 / 2, b / 2 + t_w / 2, b / 2 + t_w / 2, b / 2 + b_f2 / 2, b / 2 + b_f2 / 2, b / 2 - b_f2 / 2, b / 2 - b_f2 / 2, b / 2 - t_w / 2, b / 2 - t_w / 2, b / 2 - b_f1 / 2, b / 2 - b_f1 / 2];
    y_p = [0, 0, t_f1, t_f1, h - t_f2, h - t_f2, h, h, h - t_f2, h - t_f2, t_f1, t_f1, 0];
    
    //draw clipping mask
    svgPath = `M${x_p[0]},${y_p[0]}`
    for (let i = 1; i < x_p.length; i++) {
        segment = `L${x_p[i]},${y_p[i]}`
        svgPath = svgPath + segment
    }
    svgPath = svgPath + 'Z'
    data = { ...data, x_p: x_p, y_p: y_p, svgPath: svgPath };
    //move to geom inputs
    [y_c, width] = getWidth(data)
    data = { ...data, width: width };
    // Material properties
    E_s = 200e9; //move to geom inputs
    f_yd = data.f_yd*1e6; //move to geom inputs
    //section class 
    eps = Math.sqrt(235*1e6/f_yd)
    let SC_w = [0]
    let SC_f1 = 0
    let SC_f2 = 0
    const slenderness_w = (h-b_f1-b_f2) / t_w
    const slenderness_f1 = ((b_f1/2)-(t_w/2)) / t_f1
    const slenderness_f2 = ((b_f2/2)-(t_w/2)) / t_f2
    
    if (slenderness_w<(72*eps)) {SC_w = 1} 
    else if (slenderness_w<(83*eps)) {SC_w=2}
    else if (slenderness_w<(124*eps)) {SC_w=3}
    else {SC_w=4}
    if (slenderness_f1<(9*eps)) {SC_f1 = 1} 
    else if (slenderness_f1<(10*eps)) {SC_f1=2}
    else if (slenderness_f1<(14*eps)) {SC_f1=3}
    else {SC_f1=4}
    if (slenderness_f2<(9*eps)) {SC_f2 = 1} 
    else if (slenderness_f2<(10*eps)) {SC_f2=2}
    else if (slenderness_f2<(14*eps)) {SC_f2=3}
    else {SC_f2=4}



    eps_s = linspace(1e-4, 10e-3, 100);
    let f = (x) => x * E_s * (x <= f_yd / E_s & x >= -f_yd / E_s) + f_yd * (x > f_yd / E_s) - f_yd * (x < -f_yd / E_s);

    dy = h / y_c.length;
    y_num = y_c;
    x_num = width;
    //calculate the stress distribution, multiply by width and repeat for a lot
    //of different locations of the neutral layer. Find the solution that satisfies
    //horizontal equilibrium.
    meanStress = [];
    id = [];
    y_NA = [];
    y_tp = linspace(0.01, h, 200);
    for (let i = 0; i < eps_s.length; i++) {
        for (let j = 0; j < y_tp.length; j++) {
            strains = math.subtract(eps_s[i], math.dotMultiply(y_num, eps_s[i] / y_tp[j]));
            strains = strains.map(a => a > f_yd / E_s ? f_yd / E_s : a)
            strains = strains.map(a => a < -f_yd / E_s ? -f_yd / E_s : a)
            meanStress[j] = math.mean(math.dotMultiply(strains, x_num))
        }
        id[i] = findzero(meanStress)[0]
        y_NA[i] = y_tp[id[i]];
    }

    // now we know the NL for every strain we can easily calculate the stresses,
    // and corresponding moment. 
    const sigma_s_dist = [];
    const M = [];
    const eps_s_dist = [];
    const rinv = math.dotDivide(eps_s, y_NA);
    for (let i = 0; i < eps_s.length; i++) {
        eps_s_dist[i] = math.subtract(eps_s[i], math.dotMultiply(y_num, rinv[i]));
        sigma_s_dist[i] = math.dotMultiply(eps_s_dist[i], E_s);
        sigma_s_dist[i] = sigma_s_dist[i].map(a => a > f_yd ? f_yd : a)
        sigma_s_dist[i] = sigma_s_dist[i].map(a => a < -f_yd ? -f_yd : a)
        //M[i] = sum(math.abs(sigma_s_dist[i])*x_num*dy*abs(y_NA(i)-y_num));
        M[i] = math.sum(math.dotMultiply(math.dotMultiply(math.dotMultiply(math.abs(sigma_s_dist[i]), x_num), dy), math.abs(math.subtract(y_NA[i], y_num))))
    }
    results = { eps_s_dist: eps_s_dist, M: M, rinv: rinv, sigma_s_dist: sigma_s_dist,y_NA:y_NA };
    plot()
}



function plot() {
    ax1 = document.getElementById('myAx1')
    ax2 = document.getElementById('myAx2')
    ax3 = document.getElementById('myAx3')

    //Plot all!
    const { eps_s_dist, sigma_s_dist, M, rinv,y_NA } = results;
    const {h,b,svgPath } = data;
    val = document.getElementById('myRange').value;
    k = parseInt(val);

    let z = [];
    for (let i = 0; i < sigma_s_dist[k].length; i++) {
        z[i] = [sigma_s_dist[k][i], sigma_s_dist[k][i]]
    }

    var surf = {
            z: z,
            x: [0, b],
            y: y_c,
            type: 'contour',
            showscale: true,
            contours: {
                coloring: 'heatmap',
            },
            coloraxis: 'coloraxis',
            line: { width: 0. },
            xaxis: 'x',
            yaxis: 'y',
        }
    
    var trace21 = {
        x:[0,...sigma_s_dist[k],0],
        y:[0,...y_c,h],
        xaxis:'x2',
        yaxis:'y2',
        fill: 'toself',
        fillcolor:'rgba(0,0,0,0.2)',
        line: {width:1,color:'rgb(0,0,0)'}
    }

    var trace31 = {
        x:[0,...sigma_s_dist[k].map((x,i)=>x*width[i]/1000),0],
        y:[0,...y_c,h],
        xaxis:'x3',
        yaxis:'y3',
        fill: 'toself',
        fillcolor:'rgba(0,0,0,0.2)',
        line: {width:1,color:'rgb(0,0,0)'}
    }


    var layout = {
        paper_bgcolor: 'rgba(0,0,0,0)',
        grid: {
            rows: 1, columns: 3, pattern: 'independent',
            xgap: 0.2, ygap: 0,
            subplots: ['xy2', 'x2y2', 'x3y2'],
        },
        //coloraxis: {cmid:0,colorscale:cscale,colorbar: {x: -0.2,orientation:'h'}},
        xaxis: {
            range: [0, b],
            showgrid: false,
            showline: false,
            zeroline: false,
            
        },
        yaxis: { 
            range: [0, h], 
            scaleanchor: "x",
            showgrid: false,
            showline: false,
            zeroline: false,            
         },
        xaxis2: {
            range: [-data.f_yd*1e6*1.1, data.f_yd*1e6*1.1], 
            title: "stress (Pa)", 
            domain: [0.30, 0.5] 
        },
        yaxis2: { 
            range: [0, h], 
            showticklabels: false 
        },
        xaxis3: { 
            range: [-data.f_yd*1e6*1.1*Math.max(...width)/1000, data.f_yd*1e6*1.1*Math.max(...width)/1000], 
            domain: [0.55, 0.75],
            title: 'stress x width (N/m)'},
        yaxis3: {
            range: [0, h], showticklabels: false,
        },
        showlegend: false,
        coloraxis: {
            cmax:data.f_yd*1e6,
            cmin:-data.f_yd*1e6, 
            colorbar: {
                orientation:'v',
                thickness:15,
                x:-0.15,
                dtick:100e6,
                tickangle :45,
                ticklabelposition:"outside",
                tickcolor:'rgb(0,0,0)',
                tickfont:{
                    size:10,
                    },
                    title:{text:'(Pa)'},
                }
        }, 
        margin: {
            l: 150,
            r: 10,
            b: 50,
            t: 50,
            pad: 4
        },
        shapes: [
            {
                type: 'path',
                path: 'M-1000,-1000 L-1000,1000 L1000,1000 L1000,-1000 L-1000,-1000Z' + svgPath,
                line: {
                    color: 'rgba(0, 0, 0,1)', width: 1
                },
                fillcolor: 'rgba(255, 255, 255, 1)',
                xaxis: 'x',
                yaxis: 'y', 
            },
            {
                type: 'path',
                path: `M0,${y_NA[k]} L1,${y_NA[k]}Z`,
                line: {
                    color: 'rgba(255, 0, 0,0.5)', width: 1
                },
                xref: 'paper',
                yref: 'y1', 
            
            }],
            annotations: [{
                xref: 'paper',
                yref: 'y1',
                x: 0.9,
                xanchor: 'right',
                y: y_NA[k],
                yanchor: 'bottom',
                text: 'neutral axis',
                showarrow: false
              }]
    };
    
    Plotly.newPlot(ax1, [surf,trace21,trace31], layout,{responsive: true}  )

    var layout2 = {
        margin: {
            l: 50,
            r: 10,
            b: 50,
            t: 10,
            pad: 4
          },
        yaxis: {
            title:"Moment / kNm",
            range: [0, Math.max(...M)*1.1] },
        xaxis: {
            title:"curvature / r<sup> -1</sup>",
            range: [0, Math.max(...rinv)*1.1] },
        showlegend: false,
        paper_bgcolor:'rgba(0,0,0,0)',        
    }

    Plotly.newPlot(ax2, [{ x: rinv, y: M }],layout2)
}




calc()
document.getElementById('myRange').oninput = (x) => plot()
console.log('done')