# Fork of fTetWild

This is a fork of [`fTetWild`](https://github.com/wildmeshing/fTetWild), updated to build with modern toolchains such as **GCC ≥11** and **CMake ≥3.16**.

The goal is to ensure compatibility with current dependencies, enable CI/CD, and progressively improve stability through testing and maintenance.

---

## 📦 Dependencies Overview

| Dependency                     | Repository URL                                                                                         | GIT\_TAG             | Versioned (Release Tag)? | Notes                             |
| ------------------------------ | ------------------------------------------------------------------------------------------------------ | -------------------- | ------------------------ | --------------------------------- |
| **libigl**                     | [https://github.com/libigl/libigl.git](https://github.com/libigl/libigl.git)                           | `v2.6.0`             | ✅ Yes                    | Official release                  |
| **json**                       | [https://github.com/jdumas/json](https://github.com/jdumas/json)                                       | `0901d33...`         | ❌ No                     | No tags available                 |
| **Catch2**                     | [https://github.com/catchorg/Catch2.git](https://github.com/catchorg/Catch2.git)                       | `v3.1.1`             | ✅ Yes                    | Stable release                    |
| **CLI11**                      | [https://github.com/CLIUtils/CLI11.git](https://github.com/CLIUtils/CLI11.git)                         | `v2.5.0`             | ✅ Yes                    | Stable release                    |
| **TBB**                        | [https://github.com/uxlfoundation/oneTBB.git](https://github.com/uxlfoundation/oneTBB.git)             | `v2022.2.0-rc1`      | ⚠️ Partial               | Release candidate only            |
| **sanitizers-cmake**           | [https://github.com/arsenm/sanitizers-cmake.git](https://github.com/arsenm/sanitizers-cmake.git)       | `6947cff...`         | ❌ No                     | No version tags                   |
| **fmt**                        | [https://github.com/fmtlib/fmt.git](https://github.com/fmtlib/fmt.git)                                 | `9.0.0` *(outdated)* | ✅ Yes                    | ➜ Recommend updating to `v10.2.1` |
| **spdlog**                     | [https://github.com/gabime/spdlog.git](https://github.com/gabime/spdlog.git)                           | `v1.11.0`            | ✅ Yes                    | Stable release                    |
| **geogram**                    | [https://github.com/BrunoLevy/geogram.git](https://github.com/BrunoLevy/geogram.git)                   | `v1.9.6`             | ✅ Yes                    | Stable release                    |
| **aabbcc**                     | [https://github.com/lohedges/aabbcc.git](https://github.com/lohedges/aabbcc.git)                       | `0c85e61...`         | ❌ No                     | No tags, fixed to commit          |
| **exact envelope**             | [https://github.com/wangbolun300/fast-envelope](https://github.com/wangbolun300/fast-envelope)         | `520ee04...`         | ❌ No                     | No versioning, research code      |
| *(optional)* **WindingNumber** | [https://github.com/alecjacobson/WindingNumber.git](https://github.com/alecjacobson/WindingNumber.git) | `bde8780...`         | ❌ No                     | Unversioned, currently disabled   |

---

## 📝 Notes

* Dependencies without release tags are pinned to specific **commit hashes** to ensure reproducible builds.
* Some of these (like `aabbcc`, `fast-envelope`, `sanitizers-cmake`) have **no official versioning** and are expected to remain pinned.
* The `fmt` library should be **updated to `v10.2.1`** soon to benefit from improved C++20/C++23 support and bug fixes.

---

Would you like this to be saved as a full `README.md` file? I can output that next.
