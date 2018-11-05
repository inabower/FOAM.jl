using LinearAlgebra
import DataStructures

using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "..")
@pyimport src.toJulia as pyFoam

function readRegion(case::String)
	if !isfile(join([case,"constant","regionProperties"],"/"))
		return Dict("master"=>"None")
	end
    open(join([case,"constant","regionProperties"],"/"),"r") do f
        nameFlg = false
        bracketFlg = false
        regions = Dict()
        for line in eachline(f)
            if startswith(line, "regions")
                nameFlg = true
                continue
            end
            if nameFlg && startswith(line, "(")
                bracketFlg = true
                continue
            end
            if nameFlg && bracketFlg && startswith(line, ")")
                break
            end
            if nameFlg && bracketFlg
                #println(line)
                p0 = findfirst(isequal('('),line)
                addList = split(replace(line[p0+1:end],")"=>""))
                for add in addList
                    regions[add] = replace(line[1:p0-1]," "=>"")
                end
            end
        end
        return regions
    end
end

function timeList(case::String)
    times = Int64[]
    timeName = String[]
    for line in readdir(case)
        i = match(r"^([0-9]+)$",line)
        f = match(r"^([0-9]+\.[0-9]+)$",line)
        # println(line)
        # println(m)
        if i != nothing
            # println(i)
            if length(times) > 0
                if typeof(times[1]) == Int64
                    append!(times,[parse(Int64,i.captures[1])])
                    append!(timeName,[i.captures[1]])
                else
                    append!(times,[parse(Float64,i.captures[1])])
                    append!(timeName,[i.captures[1]])
                end
            else
                append!(times,[parse(Int64,i.captures[1])])
                append!(timeName,[i.captures[1]])
            end
        end
        if f != nothing
            # println(f)
            if length(times) > 0
                if typeof(times[1]) == Int64
                    times = convert(Array{Float64},times)
                end
            else
                times = Float64[]
            end
            append!(times,[parse(Float64,f.captures[1])])
            append!(timeName,[f.captures[1]])
        end
    end
    timelist = DataStructures.SortedDict(zip(times,timeName)) 
    return timelist
    
end

function latestTime(case::String)
    return timeList(case)[end]
end

function fieldList(case::String, time::Any, region::String)
    timeDir = convert(String,time)
    if time == "latestTime"
        timeDir = string(latestTime(case))
    end
    tmp = String[]
    masterDir = ""
    if region == nothing || region == "master"
        masterDir = joinpath(case,timeDir)
        tmp = readdir(masterDir)
    else
        masterDir = joinpath(case,timeDir,region)
        tmp = readdir(masterDir)
    end
    fieldList = String[]
    for f in tmp
		if isfile(joinpath(masterDir,f))
			append!(fieldList, [f])
		end
	end
	return fieldList
end

function readListFile(file::String, typ::String, isVector::Bool)
    open(file,"r") do f
        bracketFlg = false
        values = Any[]
        i = 1
        for line in eachline(f)
            if !bracketFlg && startswith(line, "(")
                bracketFlg = true
                continue
            end
            if bracketFlg && startswith(line, ")")
                break
            end
            if bracketFlg
                if isVector
                    p0 = findfirst(isequal('('),line)
                    addList = split(line[p0+1:end-1])
                    vec = Int64[]
                    if typ != "int"
                        vec = Float64[]
                        for add in addList
                            push!(vec, parse(Float64,add))
                        end
                    else
                        for add in addList
                            push!(vec, parse(Int64,add)+1)
                        end
                    end
                    push!(values, vec)
                else
                    if length(values) == 0
                        if typ == "int"
                            values = Int64[]
                        else
                            values = Float64[]
                        end
                    end
                    if typ != "int"
                        push!(values, parse(Float64,line))
                    else
                        push!(values, parse(Int64,line)+1)
                    end
                end
            end
        end
        return values
    end
end

function meshProperties(case::String)
    regions = readRegion(case)
    polyMeshDir = Dict()
    if length(regions) > 0
        for region in keys(regions)
			if region == "master"
				polyMeshDir[region] = join([case,"constant","polyMesh"],"/")
			else
				polyMeshDir[region] = join([case,"constant",region,"polyMesh"],"/")
			end
        end
    else
        polyMeshDir["master"] = join([case,"constant","polyMesh"],"/")
    end
    meshProp = Dict()
    meshProp["regions"] = regions
    meshProp["case"] = case
    for region in keys(polyMeshDir)
        println(region)
        masterDir = get(polyMeshDir,region,0)

        faces = readListFile(join((masterDir,"faces"),"/"), "int", true)
        points = readListFile(join((masterDir,"points"),"/"), "float", true)
        owner = readListFile(join((masterDir,"owner"),"/"), "int", false)
        neighbour = readListFile(join((masterDir,"neighbour"),"/"), "int", false)
        
        boundingBox = [[1.e10,1.e10,1.e10],[-1.e10,-1.e10,-1.e10]]
        for pi in 1:length(points)
            for x in 1:3
                boundingBox[1][x] = min(points[pi][x], boundingBox[1][x])
                boundingBox[2][x] = max(points[pi][x], boundingBox[2][x])
            end
        end
        
        nCells = max(maximum(owner),maximum(neighbour))

        fCtrs = [zeros(Float64,3) for i in 1:length(faces)]
        fAreas = [zeros(Float64,3) for i in 1:length(faces)]
        fAreaValues = zeros(Float64,length(faces))
        @threads for i in 1:length(faces)
            f = faces[i]
            sumN = zeros(Float64, 3)
            sumA = 0.0
            sumAc = zeros(Float64, 3)
            nPoints = length(f)

            if nPoints == 3
                fCtrs[i] = (1.0/3.0)*(points[f[1]] + points[f[2]] + points[f[3]])
                fAreas[i] = 0.5*cross((points[f[2]] - points[f[1]]),(points[f[3]] - points[f[1]]))
            else
                fCentre = zeros(Float64, 3)
                for fi in f
                    fCentre += points[fi] / nPoints
                end
                for pi in 1:length(f)
                    thisPoint = points[f[pi]]
                    nextPoint = points[f[(pi%nPoints)+1]]
                    c = thisPoint + nextPoint + fCentre
                    n = cross((nextPoint - thisPoint),(fCentre-thisPoint))
                    a = norm(n)
                    sumN += n
                    sumA += a
                    sumAc += a * c
                end

				if sumA < 1e-8
					fCtrs[i] = fCentre
					fAreas[i] = zeros(Float64, 3)
					fAreaValues[i] = 0.0
				else
					fCtrs[i] = sumAc/sumA/3.0
					fAreas[i] = 0.5*sumN
					fAreaValues[i] = 0.5*norm(sumN)
				end
            end
        end #end of face loop

        #Calc cell volume like src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCellCentresAndVols.C
        cellCtrs = fill(zeros(Float64,3),nCells)
        cellVols = zeros(Float64,nCells)
        cEst = fill(zeros(Float64,3),nCells)
        nCellFaces = zeros(Int8, nCells)

        @threads for facei in 1:length(owner)
            cEst[owner[facei]] += fCtrs[facei]
            nCellFaces[owner[facei]] += 1
        end
        @threads for facei in 1:length(neighbour)
            cEst[neighbour[facei]] += fCtrs[facei]
            nCellFaces[neighbour[facei]] += 1
        end
        @threads for celli in 1:nCells
            cEst[celli] /= nCellFaces[celli]
        end

        @threads for facei in 1:length(owner)
            pyr3Vol = norm(dot(fAreas[facei], (fCtrs[facei] - cEst[owner[facei]])) )
            pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[owner[facei]]
            cellCtrs[owner[facei]] += pyr3Vol*pc
            cellVols[owner[facei]] += pyr3Vol
        end
        @threads for facei in 1:length(neighbour)
            pyr3Vol = norm(dot(fAreas[facei], (fCtrs[facei] - cEst[neighbour[facei]])) )
            pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[neighbour[facei]]
            cellCtrs[neighbour[facei]] += pyr3Vol*pc
            cellVols[neighbour[facei]] += pyr3Vol
        end

        @threads for celli in 1:nCells
            if abs(cellVols[celli])<1e-16
                cellCtrs[celli] /= cellVols[celli]
            else
                cellCtrs[celli] = cEst[celli]
            end
        end
        cellVols /= 3.0
        meshProp[region] = Dict("nPoints"=>length(points), "nFaces"=>length(faces), "nCells"=>nCells, "faceCentres" => fCtrs, "faceAreas" => fAreas, "faceAreaValues" => fAreaValues, "cellCentres" => cellCtrs, "cellVolumes" => cellVols, "points"=>points,"faces"=>faces, "owner"=>owner,"neighbour"=>neighbour, "boundingBox"=>boundingBox,"boundary"=>boundary(case,string(region),fAreas))
    end
    return meshProp
end

function checkMesh(meshProp::Any)
    case = meshProp["case"]
    for region in keys(meshProp["regions"])
        println("For region ", region)
        println("  nPoints = ",meshProp[region]["nPoints"])
        println("  nFaces = ",meshProp[region]["nFaces"])
        println("  nCells = ",meshProp[region]["nCells"])
        println(" Face Area")
        println("  Min = ",minimum(meshProp[region]["faceAreaValues"]))
        println("  Max = ",maximum(meshProp[region]["faceAreaValues"]))
        println(" Volume")
        println("  Min = ",minimum(meshProp[region]["cellVolumes"]))
        println("  Max = ",maximum(meshProp[region]["cellVolumes"]))
        println("  Total = ",sum(meshProp[region]["cellVolumes"]))
        println(" Bounding Box")
        println("  Min = ",meshProp[region]["boundingBox"][1])
        println("  Max = ",meshProp[region]["boundingBox"][2])
        fields = fieldList(case, "latestTime", string(region))
        print(" Fields : \n  ")
        for f in fields
            print(f," ")
        end
        println()
        boundaryList = meshProp[region]["boundary"]
        for b in keys(boundaryList)
            println(" For Boundary ", b)
            print("    surface area : ",boundaryList[b]["magSurfaceArea"])
            println(" (vec) : ",boundaryList[b]["totalSurfaceArea"])
        end
    end
end

function runCase(case::String)
	curDir = pwd() 
	cd(case)
	run(`foamRunTutorials`)
	cd(curDir)
end

function nearestCell(p::Array{Float64}, meshProp::Any)
    ans = 0
    dist = 1.0e10
    @threads for n in length(meshProp["cellCentres"])
        if dist > norm(meshProp["cellCentres"][n] - p)
            dist = norm(meshProp["cellCentres"][n] - p)
            ans = n
        end
    end
    return ans
end

function mag(f::Array)
    ans = Float64[]
    if ndims(f) == 1
        for l in 1:length(f)
            push!(ans,abs(f[l]))
        end
    else
        for l in 1:size(f)[1]
            push!(ans,norm(f[l,:]))
        end
    end
    return ans
end

function field(case::String, name::String, time::Any, region::String)
    timeDir = convert(String,time)
    if time == "latestTime"
        timeDir = string(latestTime(case))
    end
    file = join([case,timeDir,name],"/")
    if region != nothing && region != "master"
        file = join([case,timeDir,region,name],"/")
    end
    
    return pyFoam.readDict(file)
end

function boundary(case::String, region::String, fAreas::Any)
    file = ""
    if region == nothing || region == "master"
        file = joinpath(case,"constant","polyMesh","boundary")
    else
        file = joinpath(case,"constant",region,"polyMesh","boundary")
    end
    boundaries = pyFoam.readBoundaryDict(file)
    for b in keys(boundaries)
        #println(b)
        #println(fAreas)
        sFace = boundaries[b]["startFace"]+1
        eFace = sFace + boundaries[b]["nFaces"]-1
        tA = Array{Any,1}()
        sumA = Float64[0,0,0]
        magA = 0.0
        for fi in sFace:eFace
            push!(tA,fAreas[fi])
            sumA += fAreas[fi]
            magA += norm(fAreas[fi])
        end
        boundaries[b]["surfaceArea"] = tA
        boundaries[b]["totalSurfaceArea"] = sumA
        boundaries[b]["magSurfaceArea"] = magA
    end
    return boundaries
end

function boundaryField(Uf::Any, meshProp::Any, isVector::Bool)
    Ubf = Uf["boundaryField"]
    boundaries = meshProp["boundary"]
    owner = meshProp["owner"]
    U = Dict()
    for b in keys(boundaries)
        sFace = boundaries[b]["startFace"]+1
        eFace = sFace + boundaries[b]["nFaces"]-1
        if isVector
            tU = Array{Any,1}()
        else
            tU = Float64[]
        end
        Ubff = Ubf[b]
        Utype = Ubff["type"]
        bi = 1
        #println(Utype)
        @threads for fi in sFace:eFace
            if isVector
                if in("value", keys(Ubff))
                    if length(Ubff["value"]) >= (eFace - sFace)
                        push!(tU,Ubff["value"][bi,:])
                    else
                        push!(tU,Ubff["value"])
                    end
                elseif Utype == "noSlip"
                    push!(tU,Float64[0,0,0])
                elseif Utype == "zeroGradient"
                    if length(Uf["internalField"]) >= meshProp["nCells"]
                        ownerCell = owner[fi]
                        push!(tU,Uf["internalField"][ownerCell,:])
                    else
                        push!(tU,Uf["internalField"])
                    end
                end
            else
                if in("value", keys(Ubff))
                    if length(Ubff["value"]) >= (eFace - sFace)
                        push!(tU,Ubff["value"][bi])
                    else
                        push!(tU,Ubff["value"])
                    end
                elseif Utype == "uniformTotalPressure"
                    push!(tU,Ubff["p0"][2][3])
                elseif Utype == "zeroGradient"
                    if length(Uf["internalField"]) >= meshProp["nCells"]
                        ownerCell = owner[fi]
                        push!(tU,Uf["internalField"][ownerCell])
                    else
                        push!(tU,Uf["internalField"])
                    end
                end
            end
            bi += 1
        end
        U[b] = tU
    end
    return U
end
