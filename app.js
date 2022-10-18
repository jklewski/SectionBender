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

function sectionClass(data) {
    let {h,b_f1,b_f2,t_f1,t_f2,f_yd,t_w} = data
    let SC_w = 0
    let SC_f1 = 0
    let SC_f2 = 0
    const slenderness_w = (h-b_f1-b_f2) / t_w
    const slenderness_f1 = ((b_f1/2)-(t_w/2)) / t_f1
    const slenderness_f2 = ((b_f2/2)-(t_w/2)) / t_f2  
    const eps = Math.sqrt(235*1e6/f_yd)  
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
    SC = {SC_w:SC_w,SC_f1:SC_f1,SC_f2:SC_f2}
    return SC
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
    return [y_c, width]
}

function reducedSection(data) {
    let {b,t_f1,t_f2,t_w,h,f_yd,b_f1,b_f2,E_s} = data
    //calculate gross CoG
    let A = ((b_f2-t_w)*t_f2 + (b_f1-t_w)*t_f1 + h*t_w)
    let y_tp = ((b_f1-t_w)*t_f1*t_f1/2 + 
                (b_f2-t_w)*t_f2*(h-t_f1/2) +
                (h*t_w)*h/2) / A;
    const eps = Math.sqrt(235*1e6/f_yd)    
;
    if (data.SC.SC_f2 > 3) {
    //reduce compressed flange
    let sigmaRatioF = 1;
    let k_sigmaF = (x) => (8.2/(1.05+x))*(x>=0) +
        (7.81-6.29*x+9.78*(x**2))*(x<0 && x>-1) + 
        (5.98*(1-x)**2)*(x<=-1)
    let k_sigma_in_F = k_sigmaF(sigmaRatioF)

    let lambda_p_F = ((b_f1/2-t_w/2)/t_f1) / (28.4*eps*k_sigma_in_F**0.5)
    let rhoF = 1 * (lambda_p_F<0.748) + 
        Math.min((lambda_p_F-0.188)/(lambda_p_F**2),1)*(lambda_p_F>0.748)
    let b_f2_eff = ((((b_f2-t_w)/2)*rhoF)*2+t_w)
    //calculate new TP
    A = ((b_f2_eff-t_w)*t_f2 + (b_f1-t_w)*t_f1 + h*t_w)
    y_tp = ((b_f1-t_w)*t_f1*t_f1/2 + 
                (b_f2_eff-t_w)*t_f2*(h-t_f1/2) +
                (h*t_w)*h/2) / A;
    }
    if (data.SC.SC_w > 3) {
    //reduce web
    let sigmaRatioW = (y_tp-t_f1) / (h-t_f2-(y_tp-t_f1));
    let k_sigmaW = (x) => (0.578/(0.34+x))*(x>=0) +
    (1.7-5*x+17.1*(x**2))*(x<0 && x>=-1)
    let k_sigma_in_W = k_sigmaW(sigmaRatioW)
    let lambda_p_W = ((h-t_f1-t_f2)/t_w) / (28.4*eps*k_sigma_in_W**0.5)
    let rhoW = 1 * (lambda_p_W<0.673) + 
        Math.min((lambda_p_W-0.055*(3*sigmaRatioW))/(lambda_p_W**2),1)*(lambda_p_W>=0.748)

    let h_eff_W = (h-y_tp-t_f2) * rhoW;
    let h_w_eff1 = 0.6 * h_eff_W;
    let h_w_eff2 = 0.4 * h_eff_W;
    data = {...data,h_w_eff1:h_w_eff1, h_w_eff2:h_w_eff2,b_f2_eff:b_f2_eff}    
    }
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
    // Material properties
    let E_s = 200e9; //move to geom inputs
    let f_yd = data.f_yd*1e6; //move to geom inputs
    data = { ...data, b: Math.max(data.b_f1, data.b_f2),E_s:E_s,f_yd:f_yd}
    let SC = sectionClass(data);
        
    const { h, t_f2, t_f1, t_w, b_f1, b_f2, b } = data

    data = {...data,SC:SC};
    
    //define cross section polygon
    x_p = [b / 2 - b_f1 / 2, b / 2 + b_f1 / 2, b / 2 + b_f1 / 2, b / 2 + t_w / 2, b / 2 + t_w / 2, b / 2 + b_f2 / 2, b / 2 + b_f2 / 2, b / 2 - b_f2 / 2, b / 2 - b_f2 / 2, b / 2 - t_w / 2, b / 2 - t_w / 2, b / 2 - b_f1 / 2, b / 2 - b_f1 / 2];
    y_p = [0, 0, t_f1, t_f1, h - t_f2, h - t_f2, h, h, h - t_f2, h - t_f2, t_f1, t_f1, 0];
    
    if (SC.SC_w > 3 || SC.SC_f1 > 3 || SC.SC_f2 > 3) {
        reducedSection(data);
    }
    

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
    1
    data = { ...data, width: width };
    

    



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
    const {h,b,b_f1,b_f2,t_f1,t_f2,t_w,svgPath } = data;
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


    let arrowAnnotations = [];
    let xa = [-10,b/2-b_f2/2,b/2-b_f1/2,b/2-t_w/2,b/2+t_w/2,b_f1*0.75,b_f1*0.75,b_f1*0.75,b_f1*0.75]
    let axa = [-10,b/2+b_f2/2,b/2+b_f1/2,b/2-t_w/2-20,b/2+t_w/2+20,b_f1*0.75,b_f1*0.75,b_f1*0.75,b_f1*0.75]
    let ya = [0,h+10,0-10,h/2,h/2,h,h-t_f1,0,t_f2]
    let aya = [h,h+10,0-10,h/2,h/2,h+20,h-t_f1-60,0-60,t_f2+20]
    
    for (let i=0;i<3;i++) {
        arrowAnnotations[i] = {
        showarrow:true,
        arrowhead:1,
        arrowwidth:1,
        line:{width:1},
        arrowside:"start+end",
        axref:"x",
        ayref:"y",
        yanchor: 'bottom',
        x:xa[i],
        y:ya[i],
        ay:aya[i],
        ax:axa[i]}
    }
    
    let tList = ['h','b<sub>f1<\sub>','b<sub>f2<\sub>','t<sub>f2<\sub>','t<sub>w<\sub>','t<sub>f1<\sub>'] 
    let textAnnotation = [];
    xa2 = [0,b/2,b/2,b_f1*0.75,b/2+t_w,b_f1*0.75]
    ya2 = [data.h/2,0,h,h-t_f1,h/2,t_f1]
    for (let i=0;i<6;i++) {
        textAnnotation[i] = {
            text: tList[i],
            x:xa2[i],
            y:ya2[i],
            showarrow:false,
            xanchor:'left',
            yanchor:'bottom',
        }
    }
    textAnnotation[1].xanchor = 'top'
    textAnnotation[1].yanchor = 'top'
    textAnnotation[3].yanchor = 'top'
    

    for (let i=3;i<9;i++) {
        arrowAnnotations[i] = {
            showarrow:true,
            arrowhead:1,
            arrowwidth:1,
            line:{width:1},
            arrowside:"end",
            axref:"x",
            ayref:"y",
            yanchor: 'bottom',
            x:xa[i],
            y:ya[i],
            ay:aya[i],
            ax:axa[i]}
        }
    

    pad = 30;
    var layout = {
        paper_bgcolor: 'rgba(0,0,0,0)',
        grid: {
            rows: 1, columns: 3, pattern: 'independent',
            xgap: 0.2, ygap: 0,
            subplots: ['x1y1', 'x2y1', 'x3y3'],
        },
        //coloraxis: {cmid:0,colorscale:cscale,colorbar: {x: -0.2,orientation:'h'}},
        xaxis: {
            range: [0, b],
            showgrid: false,
            showline: false,
            zeroline: false,
            
        },
        yaxis: { 
            range: [0-pad, h+pad], 
            scaleanchor: "x",
            showgrid: false,
            showline: false,
            zeroline: false,            
         },
        xaxis2: {
            range: [-data.f_yd*1.1, data.f_yd*1.1], 
            title: "stress (Pa)", 
            domain: [0.30, 0.5] 
        },
        yaxis2: { 
            range: [0-pad, h+pad], 
            showticklabels: false,
            showgrid: false,
            showline: false,
            zeroline: false,
        },
        xaxis3: { 
            range: [-data.f_yd*1.1*Math.max(...width)/1000, data.f_yd*1.1*Math.max(...width)/1000], 
            domain: [0.55, 0.75],
            title: 'stress x width (N/m)'},
        yaxis3: {
            range: [0-pad, h+pad], 
            showticklabels: false,
            showgrid: false,
            showline: false,
            zeroline: false,
        },
        showlegend: false,
        coloraxis: {
            cmax:data.f_yd,
            cmin:-data.f_yd, 
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
              },
              ...arrowAnnotations, ...textAnnotation
              
            ]
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