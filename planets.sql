-- 1. Sjekk om tabellen finnes fra før, og slett den i så fall (valgfritt, men kjekt for testing)
IF OBJECT_ID('planets', 'U') IS NOT NULL 
    DROP TABLE planets;

-- 2. Opprett tabellen med profesjonelle SQL Server-datatyper
CREATE TABLE planets (
    name NVARCHAR(50) NOT NULL, -- NVARCHAR støtter internasjonale tegn (f.eks. ÆØÅ)
    nr INT NOT NULL             -- INT brukes for heltall i SQL Server
);

-- 3. Sett inn rader (Husk: VALUES og 'enkle fnutter' rundt tekst!)
INSERT INTO planets (name, nr) 
VALUES ('Mercury', 1);

-- 4. Hent ut dataene for å bekrefte resultatet i output-vinduet
SELECT name, nr FROM planets;
