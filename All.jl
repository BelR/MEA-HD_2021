module All
using JLD, Statistics, Plots, StatsBase, HistogramThresholding, ImageContrastAdjustment, HDF5
# ocupa BRW, tamaño_cachos
function cortar_cachos( file_brw, Salida_cachos )
    variables = BRW( file_brw );
    #
    T = variables[ "Duration" ];
    SR = variables[ "SamplingRate" ];
    ξ = variables[ "Chs" ];
    Σ = variables[ "RawData" ];
    ζ = variables[ "Factor" ];
    ο = variables[ "Offset" ];
    ε, ω = tamaño_cachos( T, SR );
    # numero de espacios que ocupa el numero de segundos de duracion del experimento
    n_char = length( string( Int( floor( T/SR ) ) ) ); 
    # numero de bines
    n = size( ε, 1 );
    # numero de espacios que ocupa el numero de cachos 
    n_bins = length( string( n ) );
    #
    for i = 1:n # numero de B a cortar ( 1->4096, 1->ω )
        BIN = zeros( size( ξ, 1 ), ω ); # preallocate
        #= valores correspondientes al BIN especifico. 
        El canal 1,1 tiene el frame 1, 4097, 8193...etc =#
        β = collect( ( ε[ i, 1 ] - 1 ):1:ε[ i, 2 ] ); 
        for j = 1:ω
            # saca esos frames del machote seguido Σ, y ponlos en array en BIN
            BIN[ :, j ] = Σ[ ( β[ j ]*size( ξ, 1 ) ) + 1:( size( ξ, 1 )*β[ j + 1 ] ) ];
        end
        BIN = Float16.( ( BIN.*ζ ) .+ ο ); # conversion a volaje
        ini = lpad( string( Int(floor( ( β[ 1 ] )/SR ) ), "s" ), n_char + 1, "0" );
        eni = lpad( string( Int( floor( ( β[ end ] + 1 )/SR ) ), "s" ), n_char + 1, "0" );
        bin_time = string( ini, "-", eni );
        BINname = string( Salida_cachos, "BIN", lpad( i, n_bins, "0" ), "_", bin_time, ".jld" );   
        save( BINname, "data", BIN );
        println( lpad( i, n_bins, "0" ), " of ", n, " saved" );
    end
end
# 
function BRW( file_brw )
    brw = h5open( file_brw, "r" );
    RecVars = read( brw[ "/3BRecInfo/3BRecVars" ] );
    ζ = RecVars[ "SignalInversion" ][ ]*(RecVars["MaxVolt"][ ] - RecVars["MinVolt"][ ])/(2^RecVars["BitDepth"][ ]);
    Chs = read( brw[ "/3BRecInfo/3BMeaStreams/Raw" ] )[ "Chs" ]; 
    ξ = zeros( Int, size( Chs, 1 ), 3 );
    # Channels, [number Row Col]
    for i = 1:size( Chs, 1 ); ξ[ i, : ] = [ i, Chs[ i ].Row, Chs[ i ].Col ]; end
    variables = Dict(
        "Offset"        => RecVars[ "SignalInversion" ][ ]*RecVars[ "MinVolt" ][ ],
        "Factor"        => ζ,
        "Duration"      => Int( size( brw[ "3BData/Raw" ], 1 )/size( ξ, 1 ) ), # verdadero tiempo registrado,
        "SamplingRate"  => round( RecVars[ "SamplingRate" ][ ], digits = 1 ),
        "BRWname"       => file_brw,
        "Chs"   	    => ξ,
        "RawData"       => brw[ "3BData/Raw" ]
        );
    return variables
end
#
function tamaño_cachos( T, SR )
    # ocupa div_n_ab
    # tratando de encontrar un tamaño adecuado de cacho.
    if isinteger( SR ) # así, se tienen cachos de 1 segundo
        if isinteger( T/SR )
            n = Int.( T/SR );
        end
    elseif isinteger( T/floor( SR ) ) 
        n = Int.( T/floor( SR ) );
    elseif isinteger( T/ceil( SR ) )
        n = Int.( T/ceil( SR ) );
    else # si no, se hace un merequetengue
        div_T = div_n_ab( T );
        div_sec = div_T/SR;
        # rango de numero cachos que se quieren cortar, normalmente alto para trabajar a gusto
        hi = 4; lo = 2; # segundos
        if !isempty( div_T )
            # busca uno de los divisores de frames dentro del rango
            selected_divs = div_T[ findall( hi .>= div_sec .>= lo ) ];
        end
        if !isempty( selected_divs ) # si hubo, agarra el primero
            n = Int( T/selected_divs[ 1 ] );
        else # si no hubo, uno predeterminado ya que
            n = 60;
        end
    end
    #
    println( " Los ", ( n ), " cachos serán de ", ( ( T/SR )/n ) ," segundos. " )
    ω = Int(ceil( T/n ) ); # numero de frames finales (tamaño del cacho en frames)
    ε = zeros( Int, n, 2 ); # preallocate
    ε[ :, 1 ] = collect( 1:ω:T ); # inicio y 
    ε[ :, 2 ] = ε[ :, 1 ] .+ ω .- 1; # fin en frames de cada cacho (para cortar)
    if !isinteger( T/n )
        println(" El ultimo cacho es más chico ")
        ε[ end, 2 ] = T; # fin en frames de cada cacho (para cortar)
    end
    return ε, ω
end
#
function div_n_ab( n::Int, lo::Int64 = 1, hi::Int64 = n )                       
#=
    Divisores del numero n entre los rangos "lo" and "hi"
=#
    ρ = collect( 1:Int( floor( sqrt( n ) ) ) ); # los numeros de 1 en 1 de la raiz cuadrada 
    σ1 = findall( n.%ρ .== 0 ); # divisores de la raiz cuadrada ( residuo = 0 )
    σ2 = Int.( ( n )./( σ1 ) ); # Sacar los pares ( de 100, 2-50, 10-10, etc.)
    σ = sort( unique( vcat( σ1, σ2 ) ) ); # remover duplicados, concatenar, ordenar
    aux1 = @isdefined lo;
    aux2 = @isdefined hi;
    if aux1 && aux2
        rn = σ[ findall( hi .>= σ .>= lo ) ];
        return rn
    else
        return σ
    end
end
# funciones
# para pasar de la nomenclatura de karel a la mia
function karel2isabel( rawdata )
    dimensiones = size( rawdata );
    data4096xnframes = zeros( dimensiones[ 1 ]*dimensiones[ 2 ], dimensiones[ 3 ] );
    for i = 1:dimensiones[ 1 ]
        for j = 1:dimensiones[ 2 ]
            data4096xnframes[ ( ( i - 1 )*dimensiones[ 1 ] + j ), : ] = rawdata[ i, j, : ];
        end
    end
    return data4096xnframes
end
#
function neighborgs( center::Int, d::Int )
# obtienen la "d"-vecindad del canal "center" #
    A = reshape( 1:4096, 64, 64 );
    A = Int.( A' );
    x_c = findall( A .== center )[ ][ 2 ]
    y_c = findall( A .== center )[ ][ 1 ]
    aux = [ ( x_c - d ),( x_c + d ), ( y_c - d ), ( y_c + d ) ]
    aux[ aux .< 1 ] .= 1;
    aux[ aux .> 64 ] .= 64;
    neigh = A[ aux[ 3 ]:aux[ 4 ], aux[ 1 ]:aux[ 2 ] ];
    return neigh
end
#
function sats(data::Array, HIthr::Int, LOthr::Int)
    saturados = vcat(findall(data .>= HIthr),findall(data .<= LOthr));
    ChFrSat = zeros( Int, size( saturados, 1 ), 2 ); # preallocation
    if !isempty( saturados )
        for j = 1:size( saturados, 1 )
            ChFrSat[ j, 1 ] = saturados[ j ].I[ 1 ]; # channel
            ChFrSat[ j, 2 ] = saturados[ j ].I[ 2 ]; # Frame
        end
    end
    PSP = round( 
        ( size( ChFrSat, 1 )*100 )/( size( data, 1 )*size( data, 2 ) ), digits = 2 
        );
    println( string( " Hay ", PSP,"% de saturación." ) );
    return ChFrSat
end
#
function ESU( ch::Array, thr::Real, d::Int )
    #= Se obtienen los frames y voltajes que sobrepasan el umbral establecido. 
    Si hay eventos supraumbral a d-frames de distancia, se selecciona aquel que 
    tenga menor voltaje. =#
    OK = findall( ch .<= thr );
    if !isempty(OK)
        init = 1;
        while init == 1
            A = OK[ 1:( end - 1 ) ]; B = OK[ 2:end ];
            C = B .- A;
            D = findall( C .<= d ) .+ 1;
            if isempty( D )
                init = 0;
            else
                nook = zeros( Int, size( D, 1 ) );
                for E = 1:size( D, 1 )
                    # si el primero es menor que el segundo, quita el segundo
                    if isless( (ch[ OK[ D[ E ] ] ]), (ch[ OK[ D[ E ] - 1 ] ] ) ) 
                        nook[ E ] = D[ E ] - 1;
                    else
                        nook[ E ] = D[ E ];
                    end
                end
                if size( nook, 1 ) > 1
                    OK[ unique( nook ) ] .= 0;
                    filter!( x -> x != 0, OK );
                else
                    OK = OK[ Bool.( OK .!= OK[ nook[ 1 ] ] ) ];
                end
            end
        end
        channelC = zeros( size( OK, 1 ), 1 );
        channelC[ :, 1 ] = OK;
    else
        channelC = [ ];
    end
    return channelC
end
#
function Umbrales( W, Lo, Hi )
    W[W .== 0] .= 1
    chs = 1:size( W, 1 ); temp = log.( W ); 
    y1 = keys( countmap( temp ) ); temp = reshape( temp, 64, 64 ); 
    chs = reshape( chs, 64, 64 );
    edges, counts = build_histogram( temp, length( y1 ) );
    
    t1 = find_threshold(Otsu(), counts[1:end], edges)
    t2 = find_threshold(MinimumIntermodes(), counts[1:end], edges)
    t3 = find_threshold(Intermodes(), counts[1:end], edges)
    t4 = find_threshold(MinimumError(), counts[1:end], edges)
    t5 = find_threshold(Moments(), counts[1:end], edges)
    t6 = find_threshold(UnimodalRosin(), counts[1:end], edges)
    t7 = find_threshold(Entropy(),counts[1:end],edges)
    t8 = find_threshold(Balanced(), counts[1:end], edges)
    t9 = find_threshold(Yen(), counts[1:end], edges)
    
    l = zeros(Int,9);    l0 = 1:9;

    x1 = log.( W ); x1[ x1 .< t1 ] .= 0; 
    x1[ x1 .!= 0 ] .= 1; l[1] = length( chs[ x1 .== 1 ] )     
    x2 = log.( W ); x2[ x2 .< t2 ] .= 0; 
    x2[ x2 .!= 0 ] .= 1; l[2] = length( chs[ x2 .== 1 ] )     
    x3 = log.( W ); x3[ x3 .< t3 ] .= 0; 
    x3[ x3 .!= 0 ] .= 1; l[3] = length( chs[ x3 .== 1 ] )   
    x4 = log.( W ); x4[ x4 .< t4 ] .= 0; 
    x4[ x4 .!= 0 ] .= 1; l[4] = length( chs[ x4 .== 1 ] )
    x5 = log.( W ); x5[ x5 .< t5 ] .= 0; 
    x5[ x5 .!= 0 ] .= 1; l[5] = length( chs[ x5 .== 1 ] )     
    x6 = log.( W ); x6[ x6 .< t6 ] .= 0; 
    x6[ x6 .!= 0 ] .= 1; l[6] = length( chs[ x6 .== 1 ] )     
    x7 = log.( W ); x7[ x7 .< t7 ] .= 0; 
    x7[ x7 .!= 0 ] .= 1; l[7] = length( chs[ x7 .== 1 ] )     
    x8 = log.( W ); x8[ x8 .< t8 ] .= 0; 
    x8[ x8 .!= 0 ] .= 1; l[8] = length( chs[ x8 .== 1 ] )     
    x9 = log.( W ); x9[ x9 .< t9 ] .= 0; 
    x9[ x9 .!= 0 ] .= 1; l[9] = length( chs[ x9 .== 1 ] )     
    while isempty(l[ Hi .> l .> Lo ])
        Hi = Hi + 50;
    end
    println( string( "Al final los limites quedaron en: ", Lo,"-",Hi ) );
    metodo = Int.(l0[l .== minimum(l[ Hi .> l .> Lo ])])[1]
    if length(metodo) > 1
        metodo = Int(metodo[1]);
    end
    println(string("Seleccion: ", l[metodo]," canales..."))
    global x = zeros( size( W, 1 ) );
    if !isempty(metodo)
        if metodo == 1
            t = find_threshold(Otsu(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Otsu")
        elseif metodo == 2
            t = find_threshold(MinimumIntermodes(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("MinimumIntermodes")
        elseif metodo == 3
            t = find_threshold(Intermodes(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Intermodes")
        elseif metodo == 4
            t = find_threshold(MinimumError(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("MinimumError")
        elseif metodo == 5
            t = find_threshold(Moments(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Moments")
        elseif metodo == 6
            t = find_threshold(UnimodalRosin(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("UnimodalRosin")
        elseif metodo == 7
            t = find_threshold(Entropy(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Entropy")
        elseif metodo == 8
            t = find_threshold(Balanced(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Balanced")
        elseif metodo == 9
            t = find_threshold(Yen(), counts[1:end], edges);
            x = log.( W ); x[ x .< t ] .= 0; x[ x .!= 0 ] .= 1;
            met = ("Yen")
        end
    else
        println("No se puede tener una buena imagen con estos parametros")
        x = [ ];
        met = "error";
    end
    fig1 = heatmap(
        reshape( x, 64, 64 ),
        aspect_ratio = 1, 
        markerstrokecolor = :white,
        grid = false,
        c = :tempo,
        title = string( "Primer paso de deteccion: ", met )
        )
    println( string( "Primer paso de deteccion: ", met ) );
    return x, fig1
end
#
function AllSat( HIthr, LOthr, data )
    NOChs = [ ]; NOFrs = [ ]; # preallocation
    # Remover saturaciones positivas y negativas
    # obtener canales saturados para promediacion (Σ).
    ChFrSat = sats(data, HIthr, LOthr);
    # Aquellos saturados durante todo el bin son descartados de la lista (ChsGachos)) #
    Σ = zeros( Int, length( countmap( ChFrSat[ :, 1 ] ) ), 2 ); # preallocation
    Σ[ :, 1 ] = Int.( keys( countmap( ChFrSat[ :, 1 ] ) ) ); # que canales
    Σ[ :, 2 ] = Int.( values(countmap( ChFrSat[ :, 1 ] ) ) ); # cuantas veces
    # mas del 50% de frames saturados, se descarta el canal de la lista de promediación
    ChsGachos = Σ[ Σ[ :, 2 ] .>= Int( floor( 0.50*size( data, 2 ) ) ), 1 ]; # los gachos 
    push!( NOChs, ChsGachos ); # lista de todos los canales gachos
    #= obtener Frames saturados para promediacion (Φ).
    Aquellos saturados durante todo el bin son descartados de la lista (FrsGachos)) =#
    Φ = zeros( Int, length( countmap( ChFrSat[ :, 2 ] ) ), 2 ); # preallocation
    Φ[ :, 1 ] = Int.( keys( countmap( ChFrSat[ :, 2 ] ) ) ); # que Frames
    Φ[ :, 2 ] = Int.( values( countmap( ChFrSat[ :, 2 ] ) ) ); # cuantas veces
    # mas del 50% de canales saturados en ese frames (FrsGachos)...son gachos
    FrsGachos = Φ[ Φ[ :, 2 ] .>= Int( floor( 0.50*size( data, 1 ) ) ), 1 ]; 
    push!( NOFrs, FrsGachos ); # lista de todos los frames gachos
    #
    # Aquí se quitan los gachos de la lista de reparables
    # final list of channels
    ChFrSat = ChFrSat[ Bool.( 1 .- in.( ChFrSat[ :, 1 ], [ ChsGachos ] ) ), : ]; 
    # final list of Frames
    ChFrSat = ChFrSat[ Bool.( 1 .- in.( ChFrSat[ :, 2 ], [ FrsGachos ] ) ), : ];  
    # ahora ChFrSat contiene solo los canales y frames saturados sin los gachos.
    # Osea, los que se pueden reparar
    for l = 1:size( ChFrSat, 1 )
        Ch = ChFrSat[ l, 1 ]; # channel and
        Fr = ChFrSat[ l, 2 ]; # frame for correction
        # Vecinos del canal gacho
        NeighChs = vec( 
            reshape( neighborgs( Ch, 1 ), length( neighborgs( Ch, 1 ) ), 1 ) ); 
        filter!( e -> e≠Ch, NeighChs ) # without the center (channel of interest)
        NeighChsFr = zeros( Int, length( NeighChs ), 2 ); #preallocation
        NeighChsFr[ :, 1 ] = NeighChs; # lista de vecinos
        # cada uno en el frame a promediar
        NeighChsFr[ :, 2 ] = repeat( [ Fr ], length( NeighChs ) ); 
        # Se remueven de la lista de vecinos los canales gachos
        NeighChsFr = NeighChsFr[ 
            Bool.( 1 .- in.( NeighChsFr[ :, 1 ], [ ChsGachos ] ) ), : ]; 
        #= Para evitar reparar el (canal,frame) con sus vecinos igual de saturados se
        remueven los voltajes  superiores a los umbrales establecidos de la lista de 
        voltajes vecinos para promediacion =#
        if !isempty( NeighChsFr ) 
            NeighVoltage = data[ NeighChsFr ][ :, 1 ]; # voltejes de la vecindad
            NeighVoltage = NeighVoltage[ 
                Bool.( 1 .- ( LOthr .<=  NeighVoltage .>= HIthr ) ) 
                ]; # voltajes de la vecindad dentro de los umbrales
            if size( NeighVoltage, 1 ) >= 3 # minimo numero de vecinos para promediar
                global data[ Ch, Fr ] = mean( NeighVoltage ); 
            else
                #= Si no hay suficientes vecinos DENTRO del rango con quienes promediar,
                es mejor matarlo, creo. La otra opción sería promediar con los frames 
                inmediatos no saturados del mismo canal...pero no estoy segura. =#
                global data[ Ch, Fr ] = 0;
            end
        else
        # Si no hay suficientes vecinos at all con quienes promediar, es mejor matarlo.
            global data[ Ch, Fr ] = 0;
        end
    end
    #
    # todos los canales gachos y todos los frames gachos se vuelven 0 
    data[ ChsGachos, : ] .= 0; data[ :, FrsGachos ] .= 0; 
    #
    return data, NOChs, NOFrs
end
#
function vecindad8(punto::Array)
    Chs,ChsArray = Channels();
    puntoA4 = ( punto[ 1 ] - 1 )*64 + punto[ 2 ]
    okok = Chs[sort(filter!(x-> x != puntoA4,vec(neighborgs(puntoA4,1)))),2:3]
    okArray = [ ]
    for j = 1:size( okok, 1 )
        push!( okArray, okok[ j, : ] );
    end
    return okArray
end
#
function ComponentesSP(DatosSignados::Array)
    #Single pass method para sacar componentes disjuntos.
    lista = copy(DatosSignados)
    componentes = Set{Any}()
    while( length( lista )!=0)
        
        x = pop!(lista) #arranca el ULTIMO elemento de la lista
        listaprofundeza=Array{Int64}[]
        componentecurlab=Array{Int64}[]
        push!(listaprofundeza, x) #Pone elementos al FINAL de la lista
        push!(componentecurlab, x)    
        profundidad=0
        while ((length(listaprofundeza)!=0) && profundidad<1000)
            y=pop!(listaprofundeza)
            for v in vecindad8(y)
                if in(v, lista)
                    deleteat!(lista, indexin(Any[v], lista))
                    push!(listaprofundeza, v) 
                    profundidad+=1
                    push!(componentecurlab, v)
                end
            end
        end
        push!(componentes, componentecurlab)    
    end
    return componentes
end
#
function gruposMayorE( x, E)
    Chs, ChsArray = Channels()
    chs = Chs[ :, 1 ]; 
    maybes = chs[ Bool.( Int.( x ) ) ];
    maybesArray = Chs[ maybes, 2:3 ];
    maybesArrayArray = [ ];
    for j = 1:size( maybesArray, 1 )
        push!( maybesArrayArray, maybesArray[ j, : ] );
    end
    paragrph = ComponentesSP( maybesArrayArray );
    nparagraph = [ 0 ];
    for x in paragrph
        nparagraph = vcat( nparagraph, length( x ) );
    end
    nparagraph = nparagraph[ 2:end ]
    gruposmaybe = collect( Set( paragrph ) );
    gurposMayorE = gruposmaybe[ nparagraph .> E ];
    return gurposMayorE
end
#
function Channels()
    Chs = zeros( Int, 4096, 3 ); 
    Chs[:,1] = 1:4096; Chs[:, 2] = repeat(1:64, inner = 64); 
    Chs[:, 3] = repeat(1:64, outer = 64);
    ChsArray = Chs[ :, 2:3 ];
    return Chs, ChsArray
end
#
function grupos( E, x )
    gruposMayor = gruposMayorE( x, E )
    Chs, ChsArray = Channels( ); 
    fig2 = scatter( ( ChsArray[ :, 1 ], ChsArray[ :, 2 ] ),
            aspect_ratio = 1, 
            leg = false,
            markershape = :square,
            markercolor = :white,
            markerstrokecolor = :white,
        )
    for  j = 1:size( gruposMayor, 1 )
        for k = 1:length( gruposMayor[ j ] )
            fig2 = scatter!(
                ( gruposMayor[ j ][ k ][ 1 ], gruposMayor[ j ][ k ][ 2 ] ),
                markershape = :square,
                markercolor = length( gruposMayor[ j ] ),
                markersize = 2,
                markerstrokecolor = :white,
                title = string(
                    "Segundo paso de deteccion: grupos mayores a ", E, " elementos" ),
                titlefontsize = 8
                )
        end
    end
    return gruposMayor, fig2
end
#
function gruposA4( gruposMayor )
    coord_grupos = [ ]
    for i = 1:length( gruposMayor )
        temp = gruposMayor[ i ]
        for j = 1:length( temp )
            punto = temp[ j ];
            puntoA4 = ( punto[ 1 ] - 1 )*64 + punto[ 2 ];
            coord_grupos = vcat( coord_grupos, puntoA4 );
        end
    end
    return coord_grupos
end
#
function gruposA2( todos )
    coord_grupos = [ 0 0 ];
    for i = 1:length( todos )
        temp = todos[ i ]
        for j = 1:length( temp )
            punto = temp[ j ]
            coords = [ punto[1] punto[2] ];
            coord_grupos = vcat( coord_grupos, coords );
        end
    end
    coord_grupos = coord_grupos[2:end,:];
    return coord_grupos
end
#
function pesos( umbral, distancia, data )
    frames_all = 0; 
    w_chs = zeros( Int, size( data, 1 ) ); w_frs = zeros( Int, size( data, 2 ) );
    for C = 1:size( data, 1 )
        canal = data[ C, : ];
        eventos = Int.( ESU( canal, umbral, distancia ) );
        if !isempty( eventos )
            w_frs[ eventos ] = w_frs[ eventos ] .+ 1;
            w_chs[ C ] = length( eventos );
        else
            continue
        end
    end
    frames_all = vcat( frames_all, w_frs ) 
    return frames_all, w_chs
end
#
function find_files( path_main::String, key::String, Γ::String ) # depende de OS!!!
# buscar en el directorio Path_main, los archivos que terminen con key
    searchdir(path_main::String, key::String) = filter( x -> endswith(x, key), readdir(path_main) );
    files = searchdir( path_main, key );
    if path_main[ end ] == collect( Γ )[ 1 ]
        files = string.( path_main, files );
    else
        files = string.( path_main, Γ, files );
    end
    return files
end
#
function checkpath( workpath::String )
# checar si el folder existe, si no, hacerlo
    if isdir( workpath ) == false
        mkpath( workpath )
    end
end
#
export cortar_cachos
export BRW
export tamaño_cachos
export div_n_ab
export karel2isabel
export neighborgs
export sats
export ESU
export Umbrales
export AllSat
export vecindad8
export ComponentesSP
export gruposMayorE
export Channels
export grupos
export gruposA4
export gruposA2
export pesos
export find_files
export checkpath
end