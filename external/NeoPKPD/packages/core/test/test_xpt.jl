# Test suite for XPT (SAS Transport) Reader

using Test
using NeoPKPD

@testset "XPT Reader" begin

    @testset "XPT Types" begin
        @test isdefined(NeoPKPD, :XPTVariable)
        @test isdefined(NeoPKPD, :XPTDataset)
    end

    @testset "XPT Functions" begin
        @test isdefined(NeoPKPD, :read_xpt)
        @test isdefined(NeoPKPD, :read_cdisc_xpt)
        @test isdefined(NeoPKPD, :read_cdisc_from_xpt)
    end

    @testset "IBM to IEEE Conversion" begin
        # Test IBM 370 floating point to IEEE 754 conversion
        # This is internal but critical for XPT reading

        # Test zero
        zero_bytes = UInt8[0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]
        @test NeoPKPD.ibm_to_ieee(zero_bytes) == 0.0

        # Test missing value (all 0x2E which is '.' in EBCDIC)
        missing_bytes = fill(UInt8(0x2E), 8)
        @test isnan(NeoPKPD.ibm_to_ieee(missing_bytes))
    end

    @testset "EBCDIC to ASCII" begin
        # Test EBCDIC conversion for common characters
        # Space is 0x40 in EBCDIC
        space_bytes = UInt8[0x40]
        @test strip(NeoPKPD.ebcdic_to_ascii(space_bytes)) == ""

        # Test null termination
        null_bytes = UInt8[0xC1, 0xC2, 0x00, 0xC3]  # "AB" then null
        result = NeoPKPD.ebcdic_to_ascii(null_bytes)
        @test result == "AB"
    end

    @testset "Read String" begin
        bytes = Vector{UInt8}("HELLO   ")
        result = NeoPKPD.read_string(bytes, 1, 8)
        @test result == "HELLO"

        # Test out of bounds
        result2 = NeoPKPD.read_string(bytes, 1, 100)
        @test result2 == ""
    end

    @testset "Big Endian Integer Reading" begin
        # Test 16-bit big endian
        bytes16 = UInt8[0x01, 0x00]  # 256 in big endian
        @test NeoPKPD.read_int16_be(bytes16, 1) == 256

        bytes16_2 = UInt8[0x00, 0x05]  # 5 in big endian
        @test NeoPKPD.read_int16_be(bytes16_2, 1) == 5

        # Test 32-bit big endian
        bytes32 = UInt8[0x00, 0x00, 0x01, 0x00]  # 256 in big endian
        @test NeoPKPD.read_int32_be(bytes32, 1) == 256
    end

    @testset "XPT Variable Structure" begin
        var = XPTVariable(
            "USUBJID",
            2,  # character
            20,
            "Subject ID",
            "",
            ""
        )

        @test var.name == "USUBJID"
        @test var.ntype == 2  # character
        @test var.length == 20
        @test var.label == "Subject ID"
    end

    @testset "XPT Dataset Structure" begin
        var1 = XPTVariable("USUBJID", 2, 20, "Subject ID", "", "")
        var2 = XPTVariable("AGE", 1, 8, "Age", "", "")

        data = [
            Dict{String,Any}("USUBJID" => "001", "AGE" => 35.0),
            Dict{String,Any}("USUBJID" => "002", "AGE" => 42.0)
        ]

        dataset = XPTDataset(
            "DM",
            "Demographics",
            [var1, var2],
            data,
            2
        )

        @test dataset.name == "DM"
        @test dataset.label == "Demographics"
        @test length(dataset.variables) == 2
        @test dataset.nobs == 2
        @test length(dataset.data) == 2
    end

    @testset "XPT CDISC Integration" begin
        # Test that XPT reader functions can be called with CDISC domain readers
        @test isdefined(NeoPKPD, :read_pc_xpt)
        @test isdefined(NeoPKPD, :read_ex_xpt)
        @test isdefined(NeoPKPD, :read_dm_xpt)
    end

    @testset "Error Handling" begin
        # Test error for non-existent file
        @test_throws SystemError read_xpt("/nonexistent/path/file.xpt")
    end

end
