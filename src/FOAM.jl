using LinearAlgebra
import DataStructures

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

function checkMesh(case::String)
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
    meshProperties = Dict()
    for region in keys(polyMeshDir)
        masterDir = get(polyMeshDir,region,0)

        faces = readListFile(join((masterDir,"faces"),"/"), "int", true)
        points = readListFile(join((masterDir,"points"),"/"), "float", true)
        owner = readListFile(join((masterDir,"owner"),"/"), "int", false)
        neighbour = readListFile(join((masterDir,"neighbour"),"/"), "int", false)

        nCells = max(maximum(owner),maximum(neighbour))

        fCtrs = Any[]
        fAreas = Any[]
        fAreaValues = Float64[]
        for f in faces
            sumN = zeros(Float64, 3)
            sumA = 0.0
            sumAc = zeros(Float64, 3)
            nPoints = length(f)

            if nPoints == 3
                push!(fCtrs,(1.0/3.0)*(points[f[1]] + points[f[2]] + points[f[3]]))
                push!(fAreas,0.5*cross((points[f[2]] - points[f[1]]),(points[f[3]] - points[f[1]])))
            else
                fCentre = zeros(Float64, 3)
                for fi in f
                    fCentre += points[fi] / nPoints
                end
                for pi in 1:length(f)
                    thisPoint = points[f[pi]]
                    nextPoint = points[f[(pi%nPoints)+1]]
                    print()
                    c = thisPoint + nextPoint + fCentre
                    n = cross((nextPoint - thisPoint),(fCentre-thisPoint))
                    a = norm(n)
                    sumN += n
                    sumA += a
                    sumAc += a * c
                end

				if sumA < 1e-8
					push!(fCtrs,fCentre)
					push!(fAreas,zeros(Float64, 3))
					push!(fAreaValues,0.0)
				else
					push!(fCtrs,sumAc/sumA/3.0)
					push!(fAreas,0.5*sumN)
					push!(fAreaValues,0.5*norm(sumN))
				end
            end
        end #end of face loop

        #Calc cell volume like src/OpenFOAM/meshes/primitiveMesh/primitiveMeshCellCentresAndVols.C
        cellCtrs = fill(zeros(Float64,3),nCells)
        cellVols = zeros(Float64,nCells)
        cEst = fill(zeros(Float64,3),nCells)
        nCellFaces = zeros(Int8, nCells)

        for facei in 1:length(owner)
            cEst[owner[facei]] += fCtrs[facei]
            nCellFaces[owner[facei]] += 1
        end
        for facei in 1:length(neighbour)
            cEst[neighbour[facei]] += fCtrs[facei]
            nCellFaces[neighbour[facei]] += 1
        end
        for celli in 1:nCells
            cEst[celli] /= nCellFaces[celli]
        end

        for facei in 1:length(owner)
            pyr3Vol = norm(dot(fAreas[facei], (fCtrs[facei] - cEst[owner[facei]])) )
            pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[owner[facei]]
            cellCtrs[owner[facei]] += pyr3Vol*pc
            cellVols[owner[facei]] += pyr3Vol
        end
        for facei in 1:length(neighbour)
            pyr3Vol = norm(dot(fAreas[facei], (fCtrs[facei] - cEst[neighbour[facei]])) )
            pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst[neighbour[facei]]
            cellCtrs[neighbour[facei]] += pyr3Vol*pc
            cellVols[neighbour[facei]] += pyr3Vol
        end

        for celli in 1:nCells
            if abs(cellVols[celli])<1e-16
                cellCtrs[celli] /= cellVols[celli]
            else
                cellCtrs[celli] = cEst[celli]
            end
        end
        cellVols /= 3.0
        meshProperties[region] = Dict("nPoints"=>length(points), "nFaces"=>length(faces), "nCells"=>nCells, "faceCentres" => fCtrs, "faceAreas" => fAreas, "faceAreaValues" => fAreaValues, "cellCentres" => cellCtrs, "cellVolumes" => cellVols)
    end
    return meshProperties
end

function runCase(case::String)
	curDir = pwd() 
	cd(case)
	run(`foamRunTutorials`)
	cd(curDir)
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
    
    dimensions = Int8[]
    internalField = Float64[]
    boundaryField = Dict()
    skipL = Int64[]
    startFlg = false
    open(file,"r") do f
        lines = readlines(f)
        for l in 1:length(lines)
			if l in skipL
				continue
			end
			line = lines[l]
			#println(line)
			if startswith(line,"dimensions")
				p0 = findfirst(isequal('['),line)
				dim = split(line[p0+1:end-2])
				for d in dim
					push!(dimensions,parse(Int8,d))
				end
			end
			if startswith(line,"internalField")
				if findfirst(isequal(';'),line) != nothing
					f = match(r"^([0-9.]+)$",line)
					internalField = parse(Float64,f.captures[1])
				else
					for ll in l+3:length(lines)
						push!(skipL,ll)
						line = lines[ll]
						if startswith(line,")")
							break
						end
						
						f = match(r"^([0-9.]+)$",line)
						push!(internalField,parse(Float64,f.captures[1]))
						
						if findfirst(isequal(')'),line) != nothing
							break
						end
					end
				end
			end
			#=
			if startswith(line,"boundaryField")
				for ll in l+2:length(lines)
					if ll in skipL
						continue
					end
					push!(skipL,ll)
					line = lines[ll]
					if startswith(line,"}")
						break
					end
					if findfirst(isequal('{'),line) != nothing
						bName = lines[ll-1]
						bDict = Dict()
						for lll in ll+1:length(lines)
							push!(skipL,lll)
							if findfirst(isequal('values'),line) != nothing
								boundaryField[bName] = bDict
								break
							end
						end
					end
				end
			end
			=#
            if length(dimensions)>0 && length(internalField)>0 && length(boundaryField)>0
				break
            end
        end
    end
    return Dict("dimensions"=>dimensions,"internalField"=>internalField,"boundaryField"=>boundaryField)
end
