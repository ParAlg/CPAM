[] Remove the raw dependency on PAM and bundle it under bzlmod / cmake
[] Add tests
[] Clean up dead code; ensure the same interface as PAM
[] Check if we need std::tuple over pair (this makes it painful to swap PAM for CPAM)


