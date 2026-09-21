CaratNameByGens := function(gens, bravais)

    local grpfile, resfile, input, name, args, tori, input_args;

	if Size(gens) = 0 then
		return fail;
	fi;
	
	if Size(gens) = 1 and gens[1] = IdentityMat(Size(gens[1])) then
		if Size(gens[1]) > 1 and Size(gens[1]) < 8 then
			tori := [ "min.1.1.1.0", "min.2.1.1.0", "min.6.1.1.0", "min.15.0.1.0", "min.58.1.1.0", "min.170.0.1.0" ];
			return tori[Size(gens[1])-1];
		fi;
		return fail;
	fi;

    # get temporary file names
    grpfile := Concatenation(CaratTmpFile( "grp." ), String(IO_getpid()));
    resfile := Concatenation(CaratTmpFile( "res." ), String(IO_getpid()));

    CaratWriteMatrixFile( grpfile, gens );	

    if bravais=true then
        input_args := Concatenation( grpfile, " -Z" );
    else
        input_args := Concatenation( grpfile, "" );
    fi;
    # execute Carat program
    # args := Concatenation( grpfile, " >/dev/null 2>&1 3>&1 4>&1" );
    CaratCommand( "Name", input_args, resfile );

    # read Carat result from file, and remove temporary files
    input := InputTextFile(resfile);
    name  := ReadAllLine(input);
    CloseStream(input);
    RemoveFile( grpfile );
    RemoveFile( resfile );

    if name = fail then
        return fail;
    fi;

    name := ReplacedString(name, "\n", "");
    name := SplitString(name, " ");

    if bravais=true then
        return Concatenation(name[2], ".", name[4], ".", name[5]);
    fi;
    return Concatenation(name[2], ".", name[4], ".", name[5], ".", name[7]);

end;

CaratName := function( grp )
    if IsAffineCrystGroupOnRight(grp) then
        return CaratNameByGens( GeneratorsOfGroup( TransposedMatrixGroup(grp) ), false );
    elif IsAffineCrystGroupOnLeft(grp) then
        return CaratNameByGens( GeneratorsOfGroup(grp), false );
    elif IsIntegerMatrixGroup(grp) then
        return CaratNameByGens( GeneratorsOfGroup(grp), true );
    fi;
    return fail;
end;
