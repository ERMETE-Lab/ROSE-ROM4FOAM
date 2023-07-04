const wordList regionNames(selectRegionNames(args, runTime));
if (regionNames.size() > 1)
{
    Info << "Error: more than one region selected!" << endl;
    abort();
}
word region = regionNames[0];

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        region,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
