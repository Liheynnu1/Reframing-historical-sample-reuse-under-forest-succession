    var NDVI = image.normalizedDifference(['B8','B4']);
    
    var EVI = image.expression(
    '2.5*((NIR-Red)/(NIR+6*Red-7.5*blue+10000))',
    {
        blue: image.select('B2'),    // 0.452-0.512, BLUE
        Red: image.select('B4'),    // 0.636-0.673 μm, RED
        NIR: image.select('B8'),    //  NIR
    });
    var REPI = image.expression(
    '705+3.5*(((B4+B7)/2-B6)/(B6-B5))',
    {
        B4: image.select('B4'),    // 0.636-0.673 μm, RED
        B5: image.select('B8'),    //  NIR
        B6: image.select('B6'), 
        B7: image.select('B7'), 
    });
    var SVVI = image.expression(
    "float(sqrt((B-((B+G+R+NIR+swir1+swir2)/5))**2+(G-((B+G+R+NIR+swir1+swir2)/5))**2+(R-((B+G+R+NIR+swir1+swir2)/5))**2+(NIR-((B+G+R+NIR+swir1+swir2)/5))**2+(swir1-((B+G+R+NIR+swir1+swir2)/5))**2+(swir2-((B+G+R+NIR+swir1+swir2)/5))**2)-sqrt((NIR-((NIR+swir1+swir2)/3))**2+(swir1-((NIR+swir1+swir2)/3))**2+(swir2-((NIR+swir1+swir2)/3))**2))",
    {
                    "B": image.select("B2"),
                    "G": image.select("B3"),
                    "R": image.select("B4"),
                    "NIR": image.select("B8"),
                    "swir1": image.select("B11"),
                    "swir2": image.select("B12")
                  });
    // NDRE
    var NDRE = image.normalizedDifference(['B8', 'B7']);
    // MCARI
    var MCARI = image.expression(
    "((B5-B4)-0.2*(B5-B3)*(B5/B4))",
    {
                    "B4": image.select("B4"),
                    "B5": image.select("B5"),
                    "B3": image.select("B3")
                  });
    // MCR
    var MCR = image.expression(
    "(B8-B5-1)/(sqrt(B8-B5)+1)",
    {
                    "B5": image.select("B5"),
                    "B8": image.select("B8")
                  });
    // NDPI
    var NDPI = image.expression(
    "(B8-(0.74*B4+0.26*B11))/(B8+(0.74*B4+0.26*B11))",
    {
                    "B4": image.select("B4"),
                    "B8": image.select("B8"),
                    "B11": image.select("B11")
                  });

    // REI
    var REI = image.normalizedDifference(['B6', 'B5']);
    // IRECI
    var IRECI = image.expression(
    "(B7-B4)/(B5/B6)",
    {
                    "B4": image.select("B4"),
                    "B5": image.select("B8"),
                    "B6": image.select("B11"),
                    "B7": image.select("B7")
                  });
    var NDVIR1 = image.normalizedDifference(['B8','B5']);
    //加入波段信息
    return image.addBands(NDVI.rename("NDVI"))
              .addBands(SVVI.rename('SVVI'))
              .addBands(EVI.rename('EVI'))
              .addBands(NDRE.rename('NDRE'))
              .addBands(REPI.rename('REPI'))
              .addBands(REI.rename('REI'))
              //.addBands(MSAVI.rename('MSAVI'))
              .addBands(MCARI.rename('MCARI'))
              .addBands(NDPI.rename('NDPI'))                 
              .addBands(MCR.rename('MCR'))
              .addBands(IRECI.rename('IRECI'))
              .addBands(NDVIR1.rename('NDVIR1'));                     
};
