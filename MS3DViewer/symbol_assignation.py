
def symbol_assignation(modi):
        ## Acetylations
    if modi == 'Acetyl':
        # Orange Star
        symbol = 'star'
        color = 'rgba(255,170,12,0.8)' 
        linecolor = 'rgb(255,170,12)'
        size = 20
    elif modi=='Glutaryl':
        # Orange Triangle-Up
        symbol='triangle-up'
        color = 'rgba(255,170,12,0.8)' 
        linecolor = 'rgb(255,170,12)'
        size=10
    elif modi=='Crotonyl':
        # Orange Diamond
        symbol='diamond'
        color = 'rgba(255,170,12,0.8)' 
        linecolor = 'rgb(255,170,12)'
        size=10
    elif modi=='Butyryl':
        # Orange Trianle-Down
        symbol='triangle-down'
        color = 'rgba(255,170,12,0.8)' 
        linecolor = 'rgb(255,170,12)'
        size=10
        
    elif modi=='Propionyl':
        # Yellow Diamond
        symbol='diamond'
        color='rgba(248,230,65,0.8)'
        linecolor = 'rgb(248,230,65)'
        size=20

    elif modi=='Succinyl':
        # Light purple Triangle-Up
        symbol='triangle-up'
        color='rgba(220,140,255,0.8)'
        linecolor = 'rgb(220,140,255)'
        size=10
        
    elif modi=='Malonyl':
        # Teal Pentagon
        symbol='pentagon'
        color='rgba(110,220,255,0.8)'
        linecolor = 'rgb(110,220,255)'
        size=10
    
    elif modi =='Pyridylacetyl':
        # Bronze hexagon
        symbol = 'hexagon'
        color='rgba(191,140,70,0.8)'
        linecolor = 'rgb(191,140,70)'
        size=15
        

    ## Phosphorylations
    elif modi=='Phospho':
        # Red Circle
        symbol='circle'
        color='rgba(250,85,90,0.8)'
        linecolor = 'rgb(250,85,90)'
        size=20

    ## Methylations
    elif modi=='Methyl':
        # Green Square
        symbol='square'
        color='rgba(130,220,110,0.8)'
        linecolor = 'rgb(130,220,110)'
        size=20
    elif modi=='Dimethyl':
        # Green Diamond
        symbol='diamond'
        color='rgba(130,220,110,0.8)'
        linecolor = 'rgb(130,220,110)'
        size=15
    elif modi=='Trimethyl':
        # Green Tiangle-Up
        symbol='triangle-up'
        color='rgba(130,220,110,0.8)'
        linecolor = 'rgb(130,220,110)'
        size=15
    elif modi == 'Carbamidomethyl':
        # Light gray hexagon
        symbol = 'hexagon'
        color='rgba(200,200,200,0.8)'
        linecolor = 'rgb(200,200,200)'
        size=15

    ## Others
    
    elif modi=='Ubiquitin':
        # Purple Circle
        symbol='circle'
        color='rgba(170,110,220,0.8)'
        linecolor = 'rgb(170,110,220)'
        size=15

    elif modi=='Sumo' or 'SUMO' in modi or 'sumo' in modi or 'Sumo' in modi:
        # Blue Circle
        symbol='circle'
        color='rgba(70,175,255,0.8)'
        linecolor = 'rgb(70,175,255)'
        size=15

    elif modi=='Benzoyl':
        # Blue Star
        symbol='star'
        color='rgba(70,175,255,0.8)'
        linecolor = 'rgb(70,175,255)'
        size=15
    
    elif modi=='Citrullin':
        # Blue Diamond
        symbol='diamond'
        color='rgba(70,175,255,0.8)'
        linecolor = 'rgb(70,175,255)'
        size=15

    elif modi=='Glycosyl':
        # Pink Hexagon
        symbol='hexagon'
        color='rgba(250,120,240,0.8)'
        linecolor = 'rgb(250,120,240)'
        size=15

    elif modi=='Palmitoyl':
        # Grey Hourglass
        symbol='hourglass'
        color='rgba(177,177,177,0.8)'
        linecolor = 'rgb(177,177,177)'
        size=15

    elif modi=='Hydroxyl' or 'Hydroxy' in modi:
        # Light Blue Circle
        symbol='circle'
        color='rgba(160,210,250,0.8)'
        linecolor = 'rgb(160,210,250)'
        size=15

    elif modi=='Carboxyl':
        # Dark Green Circle
        symbol='circle'
        color='rgba(40,155,5,0.8)'
        linecolor = 'rgb(40,155,5)'
        size=15

    elif modi=='Prenyl':
        # Brown Hourglass
        symbol='hourglass'
        color='rgba(177,130,70,0.8)'
        linecolor = 'rgb(177,130,70)'
        size=15

    elif modi=='Nitrosyl':
        # Dark Blue Circle
        symbol='circle'
        color='rgba(20,80,240,0.8)'
        linecolor = 'rgb(20,80,240)'
        size=15

    elif modi=='Sulfat':
        # Teal Circle
        symbol='circle'
        color='rgba(110,220,255,0.8)'
        linecolor = 'rgb(110,220,255)'
        size=15       
    
    else:
        symbol='cross'
        color='rgba(188,188,188,0.5)'
        linecolor = 'rgb(188,188,188)'
        size=15
    
    if 'N-Term' in modi:
        symbol = 'arrow-bar-left'

    return symbol,color,linecolor,size
